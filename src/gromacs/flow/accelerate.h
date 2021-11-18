#include "gromacs/math/vectypes.h"

#ifndef MD_FLOW_FIELD_ACCELERATE
#define MD_FLOW_FIELD_ACCELERATE

struct AccelerationFlowOpts {
    AccelerationFlowOpts() {}
    AccelerationFlowOpts(const gmx::RVec rmin, const gmx::RVec rmax)
    :doLocalAcceleration{true}, rmin{rmin}, rmax{rmax} {}

    bool contains(const rvec r) const {
        for (size_t d = 0; d < DIM; ++d)
        {
            if ((r[d] < rmin[d]) || (r[d] > rmax[d]))
            {
                return false;
            }
        }

        return true;
    }

    bool doLocalAcceleration = false;
    gmx::RVec rmin, rmax;
};

#endif // MD_FLOW_FIELD_ACCELERATE