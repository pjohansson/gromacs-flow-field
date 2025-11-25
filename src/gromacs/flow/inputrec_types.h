#ifndef MD_FLOW_FIELD_IR_TYPES
#define MD_FLOW_FIELD_IR_TYPES

namespace gmx
{
namespace flow
{

//! Input options for flow field collection, read from .mdp file by grompp
//!
//! In addition to these options, `user1-grps` will contain all the simulation
//! groups that the flow field data will be sampled from.
struct FlowFieldOptions
{
    //! Whether or not to collect flow field data
    bool doFlowFieldCollection = false;

    //! Interval in steps for sampling flow field data
    int nstsample = 0;

    //! Interval in steps for averaging and writing flow field data to disk
    int nstoutput = 0;

    //! Number of flow field grid bins along x
    int nx = 0;

    //! Number of flow field grid bins along z
    int nz = 0;
};

} // namespace flow
} // namespace gmx

#endif // MD_FLOW_FIELD_IR_TYPES
