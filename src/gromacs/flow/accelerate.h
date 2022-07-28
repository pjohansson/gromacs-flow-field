#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"

#include "grid.h"
#include "inputrec_types.h"

#ifndef MD_FLOW_FIELD_ACCELERATE
#define MD_FLOW_FIELD_ACCELERATE

#ifndef NDEBUG
#define NDEBUG
#endif

namespace flow
{

//! Grid which samples the number density (1/nm^3) inside the acceleration zone
class DensityGrid : public Grid3d<float> {
public:
    DensityGrid() {}

    DensityGrid(const gmx::RVec rmin,
                const gmx::RVec rmax)
    :doDensityScaling{true},
     end{rmax}
    {
        for (size_t i = 0; i < DIM; ++i)
        {
            origin[i] = rmin[i];
            shape[i] = static_cast<int>((end[i] - origin[i]) / resolution);
            spacing[i] = (end[i] - origin[i]) / static_cast<real>(shape[i]);
        }

        _finalize();
        reset();
    }

    //! Return the y-z area of bins in the grid
    //!
    //! NOTE: Assumes that the acceleration is directed fully along the x axis
    float bin_area() noexcept
    {
        return spacing[YY] * spacing[ZZ];
    }

    //! Set the value of all bins to 0
    void reset() noexcept
    {
        for (auto& v : values)
        {
            v = 0.0;
        }
    }

    bool contains2(const rvec r) const
    {
        for (size_t i = 0; i < DIM; ++i)
        {
            if ((r[i] < origin[i]) || (r[i] > end[i]))
            {
                return false;
            }
        }

        return true;
    }

    //! Get 1 / area_density for the bin at a position
    //!
    //! Returns 0 if the position is outside the grid.
    float get(const rvec r0, const matrix box) const
    {
        rvec r;

        copy_rvec(r0, r);

        for (size_t d = 0; d < DIM; ++d)
        {
            r[d] = fmod(r[d], box[d][d]);

            while (r[d] < 0.0)
            {
                r[d] += box[d][d];
            }
        }

        if (contains2(r))
        {
            const auto value = at_pos(r);

            if (value != 0.0)
            {
                return 1.0 / value;
            }
        }

        return 0.0;
    }

    //! Whether or not to use per-density scaling
    bool doDensityScaling = false;

    //! End of grid in system coordinates (rmin + extent)
    gmx::RVec end;

    //! Target grid resolution (final per-dim resolution depends on extent)
    float resolution = 0.25;

    //! Interval in steps between updates of the local density grid
    int64_t step_update = 5000;
};


//! Options for adding acceleration to atoms inside an acceleration zone only
struct LocalAcceleration {
    //! Empty constructor (disables local acceleration)
    LocalAcceleration() {}

    //! Full constructor (enables local acceleration)
    LocalAcceleration(const LocalAccelerationOptions &opts,
                      const double                    delta_t,
                      const matrix                    box_matrix)
    :doLocalAcceleration{opts.doLocalAcceleration},
     rmin{opts.origin},
     tau{opts.tau},
     step_max_acceleration { static_cast<int64_t>(tau / delta_t) }
    {
        for (size_t d = 0; d < DIM; d++)
        {
            check_axis[d] = opts.extent[d] >= 0.0;
            box[d] = box_matrix[d][d];

            if (!check_axis[d])
            {
                rmin[d] = 0.0;
                rmax[d] = box[d];
            }
            else
            {
                rmax[d] = std::min(rmin[d] + opts.extent[d], box[d]);
            }
        }

        density_grid = DensityGrid(rmin, rmax);
    }

    /*! \brief Check if a position is inside the acceleration box.
     *
     * \param[in]   r       Position to check
     * \param[in]   box     Box size of system. Used to put r in the box before checking
     * \param[out]  result  Whether the position is inside or not
     *
     * NOTE: Always returns true if `doLocalAcceleration == false`.
     */
    bool contains(const rvec r) const
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
                p += box[d];
            }
            p = fmod(p, box[d]);

            if ((p < rmin[d]) || (p > rmax[d]))
            {
                return false;
            }
        }

        return true;
    }

    //! Whether we are doing local acceleration only or not
    bool doLocalAcceleration = false;

    //! Atom number density per area on a 3d grid
    DensityGrid density_grid;

    //! Local acceleration box origin/start
    gmx::RVec rmin;

    //! Local acceleration box end
    gmx::RVec rmax;

    //! Copy of the system box size
    //(will not change, since flow field asserts that we are not pressure scaling)
    gmx::RVec box;

    //! Input time at which full acceleration is applied (non-positive for immediate)
    real tau;

    //! Computed step at which full acceleration is applied (non-positive for immediate)
    int64_t step_max_acceleration;

    //! For each axis: whether to check the position along it
    std::array<bool, DIM> check_axis;
};

//! Calculate the current acceleration multiplier
//!
//! Uses a smooth-step function to slowly increase the acceleration
//! multiplier from 0.0 to 1.0. If step_complete <= 0, always returns 1.0.
real calc_acceleration_multiplier(const int64_t step,
                                  const int64_t step_complete);

void print_local_acceleration_info(const LocalAcceleration &opts,
                                   const gmx::MDLogger     &mdlog);

void update_local_acceleration_grid(DensityGrid            &grid,
                                    const t_commrec        *cr,
                                    const t_mdatoms        *mdatoms,
                                    const t_state          *state,
                                    const SimulationGroups *groups);

} // namespace flow

#endif // MD_FLOW_FIELD_ACCELERATE
