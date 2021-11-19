#include <array>
#include <cmath>
#include <vector>
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/topology.h"

#ifndef MD_FLOW_FIELD_ACCELERATE
#define MD_FLOW_FIELD_ACCELERATE

struct AccelerationFlowOpts {
    AccelerationFlowOpts() {}
    AccelerationFlowOpts(const t_inputrec *ir)
    :doLocalAcceleration{ir->acceleration_doLocal},
     rmin{ir->acceleration_local_origin}
    {
        // During serializing (tpxio.cpp) we always allocate DIM elements
        // for the real* arrays, thus we can safely cast to gmx::RVec
        // without worrying about out of bounds access
        const gmx::RVec extent {ir->acceleration_local_extent};

        for (size_t d = 0; d < DIM; d++)
        {
            check_axis[d] = extent[d] >= 0.0;
            rmax[d] = rmin[d] + extent[d];
        }
    }

    /*! \brief Check if a position is inside the acceleration box.
     *
     * \param[in]   r       Position to check
     * \param[in]   box     Box size of system. Used to put r in the box before checking
     * \param[out]  result  Whether the position is inside or not
     *
     * Note: Always returns true if `doLocalAcceleration == false`.
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

    real repr(const gmx::RVec& vec, const size_t d) const
    {
        return check_axis[d] ? vec[d] : -1.0;
    }

    void print_info() const
    {
        if (doLocalAcceleration)
        {
            fprintf(stderr, "\n");
            fprintf(stderr, "Local acceleration: %s\n", doLocalAcceleration ? "yes" : "no");
            fprintf(stderr, "  Origin:\t[%f, %f, %f]\n", repr(rmin, XX), repr(rmin, YY), repr(rmin, ZZ));
            fprintf(stderr, "  End:   \t[%f, %f, %f]\n", repr(rmax, XX), repr(rmax, YY), repr(rmax, ZZ));
            fprintf(stderr, "\n");
        }
    }

    //! Whether we are doing local acceleration only or not
    bool doLocalAcceleration = false;
    //! Local acceleration box origin/start
    gmx::RVec rmin;
    //! Local acceleration box end
    gmx::RVec rmax;
    //! For each axis: whether to check the position along it
    std::array<bool, DIM> check_axis;
};

#endif // MD_FLOW_FIELD_ACCELERATE