#include <array>
#include <charconv> // `from_chars`
#include <string>
#include <vector>

#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
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

void read_acceleration_opts(std::vector<t_inpfile>   &inp,
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

    printStringNoNewline(&inp, "END LOCAL ACCELERATION");
}


void check_acceleration_opts(const t_inputrec *ir,
                             WarningHandler   *wi)
{
    const auto& opts = ir->localAccelerationOptions;

    if (ir->eI != IntegrationAlgorithm::MD)
    {
        wi->addError("accelerate-local is only supported by integrator = md");
    }

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


void do_tpx_acceleration(gmx::ISerializer *serializer,
                         LocalAccelerationOptions &opts)
{
    serializer->doBool(&opts.doLocalAcceleration);
    serializer->doRvec(as_rvec_array(&opts.origin));
    serializer->doRvec(as_rvec_array(&opts.extent));
    serializer->doReal(&opts.tau);
}


void pr_acceleration(FILE*                           fp,
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


/********************
 * PRESSURE OPTIONS *
 ********************/

void read_pressure_opts(std::vector<t_inpfile>      &inp,
                        AccelerationPressureOptions &opts,
                        WarningHandler              *wi)
{
    printStringNewline(&inp, "PRESSURE (FLOW FIELD MOD.)");
    printStringNoNewline(&inp, "Apply pressure instead of group acceleration:");
    printStringNoNewline(&inp, "Treats values from accelerate as a pressure in units of kJ/mol/nm^3");
    printStringNoNewline(&inp, "and uses a local density grid to apply a constant pressure along");
    printStringNoNewline(&inp, "the chosen axis.");
    opts.doPressure = (getEnum<Boolean>(&inp, "accelerate-pressure", wi) != Boolean::No);

    printStringNoNewline(&inp, "Principal axes for 2D grid which densities are sampled on");
    opts.grid_axes = getEnum<GridAxes>(&inp, "accelerate-pressure-grid-axes", wi);

    printStringNoNewline(&inp, "Interval in steps between updating the density grid");
    opts.step_update = get_eint(&inp, "accelerate-pressure-nstupdate", 1000, wi);

    printStringNoNewline(&inp, "Resolution of grid to compute local densities for applied pressure");
    opts.density_grid_resolution = get_ereal(&inp, "accelerate-pressure-grid-resolution", 0.0, wi);
    printStringNoNewline(&inp, "Smoothing factor (sigma in a Gaussian smoothing kernel) for local density grid");
    opts.density_grid_smoothing = get_ereal(&inp, "accelerate-pressure-grid-sigma", 0.0, wi);
}


static const char* axisToString(const int axis)
{
    constexpr std::array<const char*, DIM> names = {"X", "Y", "Z"};

    if (axis < DIM)
    {
        return names.at(axis);
    }
    else
    {
        gmx_fatal(FARGS, "accelerate-pressure-axis stores invalid axis %d", axis);
    }
}


void check_pressure_opts(const t_inputrec *ir,
                         WarningHandler   *wi)
{
    const auto& opts = ir->accelerationPressureOptions;

    if (ir->eI != IntegrationAlgorithm::MD)
    {
        wi->addError("accelerate-pressure is only supported by integrator = md");
    }

    if (opts.step_update < 1)
    {
        wi->addError("accelerate-pressure-nstupdate should be >= 1");
    }
    else if (opts.step_update < 100)
    {
        const auto message = gmx::formatString(
            "accelerate-pressure-nstupdate (%d) is very low for an "
            "expensive operation, are you sure about this?",
            opts.step_update
        );
        wi->addWarning(message);
    }

    if (opts.density_grid_resolution <= 0.0)
    {
        const auto message = gmx::formatString(
            "accelerate-pressure-grid-resolution (%g) should be positive",
            opts.density_grid_resolution
        );
        wi->addError(message);
    }

    if (opts.density_grid_smoothing < 0.0)
    {
        const auto message = gmx::formatString(
            "accelerate-pressure-grid-sigma (%g) should not be negative",
            opts.density_grid_smoothing
        );
        wi->addError(message);
    }
}


void triple_check_pressure_opts(t_inputrec     *ir,
                                WarningHandler *wi)
{
    std::array<bool, DIM> acceleration_axes {false, false, false};

    // Assert that we only have acceleration along a single axis,
    // and determine which. This is used to create the 2D grid
    // to sample densities in.
    for (int n = 0; n < ir->opts.ngacc; ++n)
    {
        for (size_t i = 0; i < DIM; ++i)
        {
            if (ir->opts.acceleration[n][i] != 0.0)
            {
                acceleration_axes[i] = true;
            }
        }
    }

    int num_axes = 0;

    for (size_t i = 0; i < DIM; ++i)
    {
        if (acceleration_axes[i])
        {
            ir->accelerationPressureOptions.axis_pressure = i;
            ++num_axes;
        }
    }

    if (num_axes == 0)
    {
        wi->addError("accelerate-pressure = yes but no acceleration is set");
    }
    else if (num_axes > 1)
    {
        wi->addError(
            "accelerate-pressure = yes requires acceleration to be "
            "along a single axis for all groups"
        );
    }
}


//! Do .tpx input/output for pressure acceleration options
void do_tpx_pressure(gmx::ISerializer            *serializer,
                     AccelerationPressureOptions &opts)
{
    serializer->doBool(&opts.doPressure);
    serializer->doEnumAsInt(&opts.grid_axes);
    serializer->doInt(&opts.step_update);
    serializer->doReal(&opts.density_grid_resolution);
    serializer->doReal(&opts.density_grid_smoothing);
    serializer->doInt(&opts.axis_pressure);
}


//! Print pressure acceleration options to log
void pr_pressure(FILE*                              fp,
                 int                                indent,
                 const AccelerationPressureOptions &opts)
{
    pr_str(
        fp, indent,
        "accelerate-pressure",
        booleanValueToString(opts.doPressure)
    );
    pr_str(fp, indent, "accelerate-pressure-axis", axisToString(opts.axis_pressure));
    pr_str(fp, indent, "accelerate-pressure-grid-axes", enumValueToString(opts.grid_axes));
    pr_int(fp, indent, "accelerate-pressure-nstupdate", opts.step_update);
    pr_real(fp, indent, "accelerate-pressure-grid-resolution", opts.density_grid_resolution);
    pr_real(fp, indent, "accelerate-pressure-grid-sigma", opts.density_grid_smoothing);
}

} // namespace flow
