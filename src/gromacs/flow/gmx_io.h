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

/*************************
 * FLOW FIELD COLLECTION *
 *************************/

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


/**********************
 * LOCAL ACCELERATION *
 **********************/

//! Read acceleration options from the parameter file
void read_acceleration_opts(std::vector<t_inpfile>   &inp,
                            LocalAccelerationOptions &opts,
                            WarningHandler           *wi);

//! Verify that local acceleration options are correctly set
void check_acceleration_opts(const t_inputrec *ir,
                             WarningHandler   *wi);

//! Do .tpx input/output for local acceleration options
void do_tpx_acceleration(gmx::ISerializer *serializer,
                         LocalAccelerationOptions &opts);

//! Print local acceleration options to log
void pr_acceleration(FILE*                           fp,
                     int                             indent,
                     const LocalAccelerationOptions &opts);


/********************
 * PRESSURE OPTIONS *
 ********************/

//! Read pressure acceleration options from the parameter file
void read_pressure_opts(std::vector<t_inpfile>      &inp,
                        AccelerationPressureOptions &opts,
                        WarningHandler              *wi);

//! Verify that pressure options are correctly set
void check_pressure_opts(const t_inputrec *ir,
                         WarningHandler   *wi);

//! Verify that pressure group options are correctly set (called after indexing is done)
void triple_check_pressure_opts(t_inputrec     *ir,
                                WarningHandler *wi);

//! Do .tpx input/output for pressure acceleration options
void do_tpx_pressure(gmx::ISerializer            *serializer,
                     AccelerationPressureOptions &opts);

//! Print pressure acceleration options to log
void pr_pressure(FILE*                              fp,
                 int                                indent,
                 const AccelerationPressureOptions &opts);

} // namespace flow

#endif // MD_FLOW_FIELD_GMX_IO
