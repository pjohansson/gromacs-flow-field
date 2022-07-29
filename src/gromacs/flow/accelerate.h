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

//! Container for a density grid used for pressure acceleration
//!
//! NOTE: Grid is always the same size as the input system.
class AccelerationPressure : public Grid3d<float> {
public:
    //! Constructor using the set .mdp options and box size
    AccelerationPressure(const AccelerationPressureOptions &opts,
                         const matrix                       box);

    //! Get the scaling factor for the acceleration at the given position
    float& get_factor_at_pos(const gmx::RVec r);
    const float& get_factor_at_pos(const gmx::RVec r) const;

    //! Divide the value of all bins by their area transverse to the pressure axis
    void div_bins_by_area();

    //! Set all values in the grid to 0
    void reset();

    bool doPressure = false;

    int axis_pressure = XX;

    real sigma = 0.0;

    real target_resolution = 0.0;

    int64_t step_update = 0;

private:
    void _make_axis_1d(const size_t axis, const matrix box_matrix);

    float _bin_area() const;

    size_t _get_index_unchecked(const gmx::RVec r) const;

    size_t _get_index_along_axis(const gmx::RVec r, const size_t axis) const;
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


//! Container for acceleration modifications done by the Flow Field module
struct AccelerationFlowField {
    AccelerationFlowField(const t_inputrec *ir, const matrix box)
    :pressure{AccelerationPressure{ir->accelerationPressureOptions, box}},
     local{LocalAcceleration{ir->localAccelerationOptions, ir->delta_t, box}} {}

    //! Configuration for adding an external pressure for acceleration
    AccelerationPressure pressure;

    //! Configuration for using an acceleration zone instead of the entire system
    LocalAcceleration local;
};


//! Calculate the current acceleration multiplier
//!
//! Uses a smooth-step function to slowly increase the acceleration
//! multiplier from 0.0 to 1.0. If step_complete <= 0, always returns 1.0.
real calc_acceleration_multiplier(const int64_t step,
                                  const int64_t step_complete);

void print_local_acceleration_info(const LocalAcceleration &opts,
                                   const gmx::MDLogger     &mdlog);

void update_local_acceleration_grid(AccelerationPressure   &pressure_grid,
                                    const t_commrec        *cr,
                                    const t_mdatoms        *mdatoms,
                                    const t_state          *state,
                                    const SimulationGroups *groups);

} // namespace flow

#endif // MD_FLOW_FIELD_ACCELERATE
