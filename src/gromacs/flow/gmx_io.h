#include <vector>

#include "gmxpre.h"

#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/iserializer.h"

#include "inputrec_types.h"

#ifndef MD_FLOW_FIELD_GMX_IO
#define MD_FLOW_FIELD_GMX_IO

namespace flow
{

//! Read flow field options from the parameter file
void read_flow_field_opts(std::vector<t_inpfile> &inp,
                          FlowFieldOptions       &opts,
                          char                   *groups,
                          WarningHandler         *wi);

//! Verify that flow field options are correctly set
void check_flow_field_opts(const t_inputrec *ir,
                           WarningHandler   *wi);

//! Do .tpx input/output for flow field options
void do_tpx_flow_field(gmx::ISerializer *serializer,
                       FlowFieldOptions &opts);

//! Print flow field options to log
void pr_flow_field(FILE* fp, int indent, const FlowFieldOptions &opts);

} // namespace flow

#endif // MD_FLOW_FIELD_GMX_IO
