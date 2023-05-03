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
void read_rnemd_opts(std::vector<t_inpfile> &inp,
                     RNEMDOptions           &opts,
                     char                   *groups,
                     warninp         *wi);

//! Verify that flow field options are correctly set
void check_rnemd_opts(const t_inputrec *ir,
                      warninp   *wi,
                      const gmx::EnumerationArray<PbcType, std::string> pbcTypeNames);

//! Verify that flow field options are correctly set, after all groups have been set
void check_rnemd_groups(const t_inputrec *ir,
                        const std::vector<std::string> &group_names,
                        warninp   *wi);

//! Do .tpx input/output for flow field options
void do_tpx_rnemd(gmx::ISerializer *serializer,
                  RNEMDOptions     &opts);

//! Print flow field options to log
void pr_rnemd(FILE* fp, int indent, const RNEMDOptions &opts);

} // namespace flow

#endif // MD_FLOW_FIELD_GMX_IO
