#include <charconv> // `from_chars`
#include <string>
#include <vector>

#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/txtdump.h"

#include "gmx_io.h"

namespace flow
{

/*********************
 * UTILITY FUNCTIONS *
 *********************/

//! Read an RVec from the mdp input
static gmx::RVec read_gmx_rvec(std::vector<t_inpfile> &inp,
                               const std::string      &name,
                               const std::string      &def, // default
                               WarningHandler         *wi)
{
    gmx::RVec result = {0.0, 0.0, 0.0};

    const auto buf = setStringEntry(&inp, name, def);
    const auto values = gmx::splitString(buf);

    if (values.size() != DIM)
    {
        const auto message = gmx::formatString(
            "Expected %d values for %s, found %lu",
            DIM, name.c_str(), values.size()
        );

        wi->addError(message);
    }
    else
    {
        for (size_t i = 0; i < DIM; ++i)
        {
            try
            {
                result[i] = gmx::fromString<real>(values.at(i));
            }
            catch(gmx::GromacsException&)
            {
                const auto message = gmx::formatString(
                    "Invalid value %s in mdp file. Expected a real number.",
                    values.at(i).c_str()
                );

                wi->addError(message);
            }
        }
    }

    return result;
}


/*************************
 * FLOW FIELD COLLECTION *
 *************************/

void read_flow_field_opts(std::vector<t_inpfile> &inp,
                          FlowFieldOptions       &opts,
                          char                   *groups,
                          WarningHandler         *wi)
{
    printStringNewline(&inp, "FLOW FIELD COLLECTION");

    printStringNoNewline(&inp, "Do flow field collection: No or Yes");
    opts.doFlowFieldCollection = (getEnum<Boolean>(&inp, "flow-field", wi) != Boolean::No);

    printStringNoNewline(&inp, "This selects the subset of atoms for the flow field");
    printStringNoNewline(&inp, "collection. You can select multiple groups, in which");
    printStringNoNewline(&inp, "case the fields for all groups combined and the fields");
    printStringNoNewline(&inp, "for all individual groups are all collected and written");
    printStringNoNewline(&inp, "to disk. If no groups are selected, all atoms in the");
    printStringNoNewline(&inp, "system will be used.");
    setStringEntry(&inp, "flow-field-grps", groups, nullptr);

    printStringNoNewline(&inp, "Interval in steps between sampling flow field data");
    opts.nstsample = get_eint(&inp, "flow-nstsample", 0, wi);
    printStringNoNewline(&inp, "Interval in steps between averaging and outputting flow field data");
    opts.nstoutput = get_eint(&inp, "flow-nstoutput", 0, wi);

    printStringNoNewline(&inp, "Number of flow field grid bins along x");
    opts.nx = get_eint(&inp, "flow-nx", 0, wi);
    printStringNoNewline(&inp, "Number of flow field grid bins along z");
    opts.nz = get_eint(&inp, "flow-nz", 0, wi);
}


void check_flow_field_opts(const t_inputrec *ir,
                           WarningHandler   *wi)
{
    const auto& opts = ir->flowFieldOptions;

    if (opts.nx < 1)
    {
        wi->addError("flow-nx should be >= 1");
    }
    if (opts.nz < 1)
    {
        wi->addError("flow-nz should be >= 1");
    }
    if (opts.nstsample < 1)
    {
        wi->addError("flow-nstsample should be >= 1");
    }
    if (opts.nstoutput < 1)
    {
        wi->addError("flow-nstoutput should be >= 1");
    }

    if (opts.nstoutput % opts.nstsample != 0)
    {
        const std::string message =
            gmx::formatString(
                "flow-nstoutput (%d) should be "
                "a multiple of flow-nstsample (%d)",
                opts.nstoutput, opts.nstsample
            );

        wi->addError(message);
    }

    if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
    {
        const std::string message =
            gmx::formatString(
                "Pressure scaling and flow field collection were both turned "
                "on, but flow field collection requires a fixed system box!"
            );

        wi->addError(message);
    }
}


void do_tpx_flow_field(gmx::ISerializer *serializer,
                       FlowFieldOptions &opts)
{
    serializer->doBool(&opts.doFlowFieldCollection);
    serializer->doInt(&opts.nstsample);
    serializer->doInt(&opts.nstoutput);
    serializer->doInt(&opts.nx);
    serializer->doInt(&opts.nz);
}


void pr_flow_field(FILE* fp, int indent, const FlowFieldOptions &opts)
{
    pr_str(
        fp, indent,
        "flow-field",
        booleanValueToString(opts.doFlowFieldCollection)
    );
    pr_int(fp, indent, "flow-nstsample", opts.nstsample);
    pr_int(fp, indent, "flow-nstoutput", opts.nstoutput);
    pr_int(fp, indent, "flow-nx", opts.nx);
    pr_int(fp, indent, "flow-nz", opts.nz);
}


/**********************
 * LOCAL ACCELERATION *
 **********************/

void read_local_acceleration_opts(std::vector<t_inpfile>   &inp,
                                  LocalAccelerationOptions &opts,
                                  WarningHandler           *wi)
{
    printStringNewline(&inp, "LOCAL ACCELERATION (FLOW FIELD MOD.)");
    printStringNoNewline(&inp, "Adds group acceleration (defined above) to atoms only when ");
    printStringNoNewline(&inp, "they are in a specified acceleration zone.");
    opts.doLocalAcceleration = (getEnum<Boolean>(&inp, "accelerate-local", wi) != Boolean::No);

    printStringNoNewline(&inp, "Area begins at an origin and is of a system absolute size (extent)");
    printStringNoNewline(&inp, "Negative extent along any dimension means use entire length");

    opts.origin = read_gmx_rvec(inp, "accelerate-local-origin", "0 0 0", wi);
    opts.extent = read_gmx_rvec(inp, "accelerate-local-extent", "0 0 0", wi);

    printStringNoNewline(&inp, "For tau positive: increase acceleration from 0 at t=0 to full at t=tau");
    opts.tau = get_ereal(&inp, "accelerate-tau", 0.0, wi);

    printStringNewline(&inp, "END LOCAL ACCELERATION");
}


void check_local_acceleration_opts(const t_inputrec *ir,
                                   WarningHandler   *wi)
{
    const auto& opts = ir->localAccelerationOptions;

    for (size_t i = 0; i < DIM; ++i)
    {
        if (opts.origin[i] < 0.0)
        {
            const auto message = gmx::formatString(
                "accelerate-local-origin (%g, %g, %g) does not lie "
                "within the box",
                opts.origin[XX], opts.origin[YY], opts.origin[ZZ]
            );

            wi->addWarning(message);
            break;
        }
    }

    for (size_t i = 0; i < DIM; ++i)
    {
        if (opts.extent[i] == 0.0)
        {
            const auto message = gmx::formatString(
                "accelerate-local-extent (%g, %g, %g) is 0 along an axis",
                opts.extent[XX], opts.extent[YY], opts.extent[ZZ]
            );

            wi->addError(message);
            break;
        }
    }

    if (opts.tau < 0.0)
    {
        const auto message = gmx::formatString(
            "accelerate-tau (%g) should be >= 0",
            opts.tau
        );

        wi->addError(message);
    }
}


void do_tpx_local_acceleration(gmx::ISerializer *serializer,
                               LocalAccelerationOptions &opts)
{
    serializer->doBool(&opts.doLocalAcceleration);
    serializer->doRvec(as_rvec_array(&opts.origin));
    serializer->doRvec(as_rvec_array(&opts.extent));
    serializer->doReal(&opts.tau);
}


void pr_local_acceleration(FILE*                           fp,
                           int                             indent,
                           const LocalAccelerationOptions &opts)
{
    pr_str(
        fp, indent,
        "accelerate-local",
        booleanValueToString(opts.doLocalAcceleration)
    );
    pr_rvec(fp, indent + 2, "origin", opts.origin, DIM, true);
    pr_rvec(fp, indent + 2, "extent", opts.extent, DIM, true);
    pr_real(fp, indent, "accelerate-tau", opts.tau);
}

} // namespace flow
