#include <array>
#include <cmath>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"

#include "inputrec_types.h"

#ifndef MD_FLOW_FIELD_ACCELERATE
#define MD_FLOW_FIELD_ACCELERATE

namespace flow
{

//! Options for adding acceleration to atoms inside an acceleration zone only
struct LocalAcceleration {
    //! Empty constructor (disables local acceleration)
    LocalAcceleration() {}

    //! Full constructor (enables local acceleration)
    LocalAcceleration(const LocalAccelerationOptions &opts,
                         const double                    delta_t)
    :doLocalAcceleration { opts.doLocalAcceleration },
     rmin { opts.origin },
     tau { opts.tau },
     step_full_acceleration { static_cast<int64_t>(tau / delta_t) }
    {
        for (size_t d = 0; d < DIM; d++)
        {
            check_axis[d] = opts.extent[d] >= 0.0;
            rmax[d] = rmin[d] + opts.extent[d];
        }
    }

    /*! \brief Check if a position is inside the acceleration box.
     *
     * \param[in]   r       Position to check
     * \param[in]   box     Box size of system. Used to put r in the box before checking
     * \param[out]  result  Whether the position is inside or not
     *
     * NOTE: Always returns true if `doLocalAcceleration == false`.
     */
    bool contains(const rvec r, const matrix box) const
    {
        if (!doLocalAcceleration)
        {
            return true;
        }

        for (size_t d = 0; d < DIM; d++)
        {
            if (!check_axis[d])
            {
                continue;
            }

            auto p = r[d];
            while (p < 0.0)
            {
                p += box[d][d];
            }
            p = fmod(p, box[d][d]);

            if ((p < rmin[d]) || (p > rmax[d]))
            {
                return false;
            }
        }

        return true;
    }

    //! For non-zero tau, calculate the current acceleration multiplier
    real calc_acceleration_multiplier(const int64_t step) const
    {
        if ((step >= step_full_acceleration) || (step_full_acceleration <= 0))
        {
            return 1.0;
        }
        else
        {
            // Smoothstep function for range [0.0, 1.0) -> [0.0, 1.0)
            const real x = static_cast<real>(step) / static_cast<real>(step_full_acceleration);
            return 6.0 * powf(x, 5.0) - 15.0 * powf(x, 4.0) + 10 * powf(x, 3.0);
        }
    }

    //! Whether we are doing local acceleration only or not
    bool doLocalAcceleration = false;

    //! Local acceleration box origin/start
    gmx::RVec rmin;

    //! Local acceleration box end
    gmx::RVec rmax;

    //! Input time at which full acceleration is applied (non-positive for immediate)
    real tau;

    //! Computed step at which full acceleration is applied (non-positive for immediate)
    int64_t step_full_acceleration;

    //! For each axis: whether to check the position along it
    std::array<bool, DIM> check_axis;
};

void print_local_acceleration_info(const LocalAcceleration &opts,
                                   const gmx::MDLogger     &mdlog);

} // namespace flow

#endif // MD_FLOW_FIELD_ACCELERATE
