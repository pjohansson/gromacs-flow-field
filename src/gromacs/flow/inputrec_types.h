#ifndef GMX_FLOW_INPUTREC_TYPES_H
#define GMX_FLOW_INPUTREC_TYPES_H

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
    int nstSample = 0;

    //! Interval in steps for averaging and writing flow field data to disk
    int nstOutput = 0;

    //! Number of flow field grid bins along x
    int numBinsX = 0;

    //! Number of flow field grid bins along z
    int numBinsZ = 0;
};

} // namespace flow
} // namespace gmx

#endif // GMX_FLOW_INPUTREC_TYPES_H
