#ifndef MD_FLOW_FIELD_IR_TYPES
#define MD_FLOW_FIELD_IR_TYPES

#include "gromacs/math/vectypes.h"

namespace flow
{

//! Input options for flow field collection, read from .mdp file by grompp
//!
//! In addition to these options, `user1-grps` will contain all the simulation
//! groups that the flow field data will be sampled from.
struct FlowFieldOptions {
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


struct LocalAccelerationOptions {
    bool doLocalAcceleration = false;

    gmx::RVec origin = {0.0, 0.0, 0.0};

    gmx::RVec extent = {0.0, 0.0, 0.0};

    real tau = 0.0;
};


enum class GridAxes : int {
    XY,
    XZ,
    YZ,
    Count,
    Default = XZ
};


struct AccelerationPressureOptions {
    bool doPressure = false;

    int axis_pressure = XX;

    GridAxes grid_axes;

    int step_update = 1000;

    real density_grid_resolution = 0.1;

    real density_grid_smoothing = 0.2;
};

} // namespace flow

// Must be outside of the flow namespace because it's dynamically accessed
// by other Gromacs functions (through functions in readinp.h)
const char* enumValueToString(flow::GridAxes enumValue);

#endif // MD_FLOW_FIELD_IR_TYPES
