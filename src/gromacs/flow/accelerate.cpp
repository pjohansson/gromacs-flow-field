#include <cmath>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/vec.h"

#include "accelerate.h"
#include "utils.h"

namespace flow
{

real calc_acceleration_multiplier(const int64_t step,
                                  const int64_t step_complete)
{
    if ((step >= step_complete) || (step_complete <= 0))
    {
        return 1.0;
    }
    else
    {
        const auto x = static_cast<real>(step) / static_cast<real>(step_complete);

        // Smoothstep function for range [0.0, 1.0) -> [0.0, 1.0)
        return 6.0 * powf(x, 5.0) - 15.0 * powf(x, 4.0) + 10 * powf(x, 3.0);
    }
}


//! Sum the number of atoms on grid and send to all ranks
static void mpi_collect_grid(DensityGrid &grid, const t_commrec *cr)
{
    if (PAR(cr))
    {
        MPI_Allreduce(
            MPI_IN_PLACE,
            grid.values.data(),
            grid.values.size(),
            MPI_FLOAT,
            MPI_SUM,
            cr->mpi_comm_mygroup
        );
    }
}


static void collect_grid_data(DensityGrid            &grid,
                              const t_commrec        *cr,
                              const t_mdatoms        *mdatoms,
                              const t_state          *state,
                              const SimulationGroups *groups)
{
    rvec r;

    const auto num_groups =
        get_num_groups(groups, SimulationAtomGroupType::Acceleration);

    for (size_t i = 0; i < static_cast<size_t>(mdatoms->homenr); ++i)
    {
        // Check for match to the input group using the global atom index,
        // since groups contain these indices instead of MPI rank local indices
        const auto global_atom_index = haveDDAtomOrdering(*cr)
            ? cr->dd->globalAtomIndices[i]
            : static_cast<int>(i);

        const auto group_index = getGroupType(
            *groups, SimulationAtomGroupType::User1, global_atom_index
        );

        if (group_index < static_cast<int>(num_groups))
        {
            copy_rvec(state->x[i], r);

            for (size_t d = 0; d < DIM; ++d)
            {
                r[d] = fmod(r[d], state->box[d][d]);

                while (r[d] < 0.0)
                {
                    r[d] += state->box[d][d];
                }
            }

            if (grid.contains(r))
            {
                grid.at_pos(r) += 1.0;
            }
        }
    }
}


class GaussianKernel {
public:
    GaussianKernel(const gmx::RVec spacing,
           const float     sigma,
           const float     cutoff)
    :half_distance{
        static_cast<int>(cutoff / spacing[XX]),
        static_cast<int>(cutoff / spacing[YY]),
        static_cast<int>(cutoff / spacing[ZZ])
    }
    {
        const auto row = std::vector<float>(2 * half_distance[ZZ] + 1, 0.0);
        const auto layer = std::vector<std::vector<float>>(2 * half_distance[YY] + 1, row);
        weights = std::vector<std::vector<std::vector<float>>>(2 * half_distance[XX] + 1, layer);

        for (int ix = -half_distance[XX]; ix <= half_distance[XX]; ++ix)
        {
            for (int iy = -half_distance[YY]; iy <= half_distance[YY]; ++iy)
            {
                for (int iz = -half_distance[ZZ]; iz <= half_distance[ZZ]; ++iz)
                {
                    const auto dx = static_cast<float>(ix) * spacing[XX];
                    const auto dy = static_cast<float>(iy) * spacing[YY];
                    const auto dz = static_cast<float>(iz) * spacing[ZZ];

                    const auto radius = std::sqrt(dx * dx + dy * dy + dz * dz);

                    if (radius <= cutoff)
                    {
                        const float weight = std::exp(-radius / (2.0 * sigma * sigma));
                        at(ix, iy, iz) = weight;
                    }
                }
            }
        }
    }

    float& at(const int ix, const int iy, const int iz)
    {
        const auto ix_adjusted = ix + half_distance[XX];
        const auto iy_adjusted = iy + half_distance[YY];
        const auto iz_adjusted = iz + half_distance[ZZ];

        return weights.at(ix_adjusted).at(iy_adjusted).at(iz_adjusted);
    }

    const float& at(const int ix, const int iy, const int iz) const
    {
        const auto ix_adjusted = ix + half_distance[XX];
        const auto iy_adjusted = iy + half_distance[YY];
        const auto iz_adjusted = iz + half_distance[ZZ];

        return weights.at(ix_adjusted).at(iy_adjusted).at(iz_adjusted);
    }

    void add_weights(const DensityGrid  &grid,
                     std::vector<float> &result,
                     std::vector<float> &weights,
                     const int           ix0,
                     const int           iy0,
                     const int           iz0)
    {
        const auto i0 = iz0 + iy0 * grid.shape[ZZ] + ix0 * grid.shape[YY] * grid.shape[ZZ];

        for (int dx = -xs(); dx <= xs(); ++dx)
        {
            const auto ix = ix0 + dx;

            if ((ix < 0) || (ix >= grid.shape[XX]))
            {
                continue;
            }

            for (int dy = -ys(); dy <= ys(); ++dy)
            {
                const auto iy = iy0 + dy;

                if ((iy < 0) || (iy >= grid.shape[YY]))
                {
                    continue;
                }

                for (int dz = -zs(); dz <= zs(); ++dz)
                {
                    const auto iz = iz0 + dz;
                    const auto i = iz + iy * grid.shape[ZZ] + ix * grid.shape[YY] * grid.shape[ZZ];

                    if ((iz < 0) || (iz >= grid.shape[ZZ]))
                    {
                        continue;
                    }

                    const auto weight = at(dx, dy, dz);

                    result.at(i0) += weight * grid.values.at(i);
                    weights.at(i0) += weight;
                }
            }
        }
    }

    int xs() const { return half_distance[XX]; }
    int ys() const { return half_distance[YY]; }
    int zs() const { return half_distance[ZZ]; }

private:
    gmx::IVec half_distance;

    std::vector<std::vector<std::vector<float>>> weights;
};


//! Use a 3d Guassian kernel to smooth the density grid
//!
//! This is likely a very expensive operation, but we shouldn't be
//! updating the local density grid very often which means that it
//! should be negligible compared to the force calculation.
static void smooth_gaussian_kernel(DensityGrid &grid,
                                   const real   sigma)
{
    if (sigma <= 0.0)
    {
        return;
    }

    // 3 sigma ~ 97 percent of all possible weights, by far sufficient
    // for this silly little smoothing kernel
    const auto cutoff = 3.0 * sigma;

    std::vector<float> result(grid.values.size(), 0.0);
    std::vector<float> weights(grid.values.size(), 0.0);

    GaussianKernel kernel(grid.spacing, sigma, cutoff);

    for (int ix = 0; ix < grid.shape[XX]; ++ix)
    {
        for (int iy = 0; iy < grid.shape[YY]; ++iy)
        {
            for (int iz = 0; iz < grid.shape[ZZ]; ++iz)
            {
                kernel.add_weights(grid, result, weights, ix, iy, iz);
            }
        }
    }

    for (size_t i = 0; i < result.size(); ++i)
    {
        if (weights.at(i) != 0.0)
        {
            grid.values.at(i) = result.at(i) / weights.at(i);
        }
    }
}


void update_local_acceleration_grid(DensityGrid            &grid,
                                    const t_commrec        *cr,
                                    const t_mdatoms        *mdatoms,
                                    const t_state          *state,
                                    const SimulationGroups *groups)
{
    grid.reset();
    collect_grid_data(grid, cr, mdatoms, state, groups);
    mpi_collect_grid(grid, cr);

    for (auto& v : grid.values)
    {
        v /= grid.bin_area();
    }

    smooth_gaussian_kernel(grid, 0.25);
}

void print_local_acceleration_info(const flow::LocalAcceleration &opts,
                                   const gmx::MDLogger           &mdlog)
{
    // Log to warning level, which prints both to md.log and stdout
    // (info level only writes to md.log)

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendText("**********************************\n")
        .appendText("* LOCAL ACCELERATION INFORMATION *\n")
        .appendText("**********************************");

    const auto& rmin = opts.rmin;
    const auto& rmax = opts.rmax;
    const auto& axis = opts.check_axis;

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendTextFormatted(
            "Local acceleration is enabled: %s\n",
            opts.doLocalAcceleration ? "yes" : "no"
        );

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendTextFormatted(
            "Origin:  [%g, %g, %g]\n",
            axis[XX] ? rmin[XX] : -1,
            axis[YY] ? rmin[YY] : -1,
            axis[ZZ] ? rmin[ZZ] : -1
        )
        .appendTextFormatted(
            "End:     [%g, %g, %g]\n",
            axis[XX] ? rmax[XX] : -1,
            axis[YY] ? rmax[YY] : -1,
            axis[ZZ] ? rmax[ZZ] : -1
        );

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendTextFormatted(
            "Activation time: %g",
            opts.tau
        );

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendText("**************************************\n")
        .appendText("* END LOCAL ACCELERATION INFORMATION *\n")
        .appendText("**************************************");
}

} // namespace flow
