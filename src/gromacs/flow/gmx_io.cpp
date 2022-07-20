#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/txtdump.h"

#include "gmx_io.h"

namespace flow
{

void read_flow_field_opts(std::vector<t_inpfile> *inp, 
                          FlowFieldOptions       &opts, 
                          WarningHandler         *wi)
{
    printStringNewline(inp, "FLOW FIELD COLLECTION");

    printStringNoNewline(inp, "Do flow field collection: No or Yes");
    opts.doFlowFieldCollection = (getEnum<Boolean>(inp, "flow-field", wi) != Boolean::No);

    printStringNoNewline(inp, "Interval in steps between sampling flow field data");
    opts.nstsample = get_eint(inp, "flow-nstsample", 0, wi);
    printStringNoNewline(inp, "Interval in steps between averaging and outputting flow field data");
    opts.nstoutput = get_eint(inp, "flow-nstoutput", 0, wi);

    printStringNoNewline(inp, "Number of flow field grid bins along x");
    opts.nx = get_eint(inp, "flow-nx", 0, wi);
    printStringNoNewline(inp, "Number of flow field grid bins along z");
    opts.nz = get_eint(inp, "flow-nz", 0, wi);
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

void pr_flow_field(FILE* fp, int indent, const flow::FlowFieldOptions &opts)
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

} // namespace flow
