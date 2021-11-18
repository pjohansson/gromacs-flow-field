#include <cmath>
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/topology.h"

#ifndef MD_FLOW_FIELD_ACCELERATE
#define MD_FLOW_FIELD_ACCELERATE

struct AccelerationFlowOpts {
    AccelerationFlowOpts() {}
    AccelerationFlowOpts(const gmx::RVec rmin, const gmx::RVec rmax)
    :doLocalAcceleration{true}, rmin{rmin}, rmax{rmax} {}

    /*! \brief Check if a position is inside the acceleration box.
     * 
     * Note: Always returns true if `doLocalAcceleration == false`.
     *
     * \param[in]   r       Position to check
     * \param[in]   box     Box size of system. Used to put r in the box before checking
     * \param[out]  result  Whether the position is inside or not
     */
    bool contains(const rvec r, const matrix box) const 
    {
        if (!doLocalAcceleration)
        {
            return true;
        }

        for (size_t d = 0; d < DIM; ++d)
        {
            auto p = fmod(r[d], box[d][d]);
            while (p < 0.0)
            {
                p += box[d][d];
            }

            if ((p < rmin[d]) || (p > rmax[d]))
            {
                return false;
            }
        }

        return true;
    }

    //! Whether we are doing local acceleration only or not
    bool doLocalAcceleration = false;
    //! Local acceleration box origin/start
    gmx::RVec rmin;
    //! Local acceleration box end
    gmx::RVec rmax;
};

#endif // MD_FLOW_FIELD_ACCELERATE