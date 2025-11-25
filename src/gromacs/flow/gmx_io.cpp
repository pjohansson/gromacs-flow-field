#include "gmx_io.h"

#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/txtdump.h"

#include "inputrec_types.h"

namespace gmx
{
namespace flow
{

void read_flow_field_opts(std::vector<t_inpfile>& inp, FlowFieldOptions& opts, char* groups, WarningHandler* wi)
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


void check_flow_field_opts(const t_inputrec* ir, WarningHandler* wi)
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
        const std::string message = gmx::formatString(
                "flow-nstoutput (%d) should be "
                "a multiple of flow-nstsample (%d)",
                opts.nstoutput,
                opts.nstsample);

        wi->addError(message);
    }

    if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
    {
        const std::string message = gmx::formatString(
                "Pressure scaling and flow field collection were both turned "
                "on, but flow field collection requires a fixed system box!");

        wi->addError(message);
    }
}


void do_tpx_flow_field(gmx::ISerializer* serializer, FlowFieldOptions& opts)
{
    serializer->doBool(&opts.doFlowFieldCollection);
    serializer->doInt(&opts.nstsample);
    serializer->doInt(&opts.nstoutput);
    serializer->doInt(&opts.nx);
    serializer->doInt(&opts.nz);
}


void pr_flow_field(FILE* fp, int indent, const FlowFieldOptions& opts)
{
    pr_str(fp, indent, "flow-field", booleanValueToString(opts.doFlowFieldCollection));
    pr_int(fp, indent, "flow-nstsample", opts.nstsample);
    pr_int(fp, indent, "flow-nstoutput", opts.nstoutput);
    pr_int(fp, indent, "flow-nx", opts.nx);
    pr_int(fp, indent, "flow-nz", opts.nz);
}

} // namespace flow
} // namespace gmx
