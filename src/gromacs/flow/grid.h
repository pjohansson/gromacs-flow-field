#ifndef MD_FLOW_FIELD_GRID
#define MD_FLOW_FIELD_GRID

#include <optional>
#include <stdexcept>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/cstringutil.h"

namespace gmx
{
namespace flow
{

template<typename T = double>
class Grid3d
{
public:
    //! Empty constructor
    Grid3d() {}

    //! Default constructor of grid, does not set values
    Grid3d(const IVec shape, const RVec spacing, const std::optional<RVec> opt_origin) :
        shape_{ shape }, spacing_{ spacing }
    {
        if (opt_origin)
        {
            origin_ = opt_origin.value();
        }

        _finalize();
    }

    //! Assign `value` to all cells in grid
    void assign(const T value)
    {
        for (auto& v : values_)
        {
            v = value;
        }
    }

    //! Grid cell reference accessors
    //!
    //! Throws `std::out_of_range` if 3D position is not within grid.
    const T& at(int ix, int iy, int iz) const;
    T&       at(int ix, int iy, int iz);

    //! Grid cell reference accessors from particle positions
    //!
    //! Cell access saturates at the grid edges: positions outside
    //! access the cell at the closest edge.
    const T& at_pos(const rvec r) const;
    T&       at_pos(const rvec r);

    //! Grid cell reference accessors from particle positions, pbc correction
    //!
    //! Cell access saturates at the grid edges: positions outside
    //! access the cell at the closest edge.
    const T& at_pos_pbc(const rvec r, const matrix box) const;
    T&       at_pos_pbc(const rvec r, const matrix box);

    //! Check whether a position is contained within the grid
    bool contains(const rvec r) const
    {
        return (r[XX] >= origin_[XX] && r[XX] <= origin_[XX] + box_[XX] && r[YY] >= origin_[YY]
                && r[YY] <= origin_[YY] + box_[YY] && r[ZZ] >= origin_[ZZ]
                && r[ZZ] <= origin_[ZZ] + box_[ZZ]);
    }

    //! Return the bin volume
    real bin_volume() const noexcept { return spacing_[XX] * spacing_[YY] * spacing_[ZZ]; }

    //! Return the grid shape.
    const IVec& shape() const { return shape_; }

    //! Return the grid spacing.
    const RVec& spacing() const { return spacing_; }

    const RVec& origin() const { return origin_; }

    const RVec& invSpacing() const { return invSpacing_; }

    std::vector<T>& values() { return values_; }

    const std::vector<T>& values() const { return values_; }

protected:
    //! Set-up the data after initializing primary variables
    //!
    //! To be called by a constructor which sets `shape` and `spacing`.
    void _finalize()
    {
        if ((shape_[XX] < 1) || (shape_[YY] < 1) || (shape_[ZZ] < 1))
        {
            char buf[STRLEN];
            snprintf(buf,
                     STRLEN,
                     "Grid3d::Grid3d: shape (%d, %d, %d) "
                     "must be positive along all directions",
                     shape_[XX],
                     shape_[YY],
                     shape_[ZZ]);

            throw std::invalid_argument(buf);
        }

        for (size_t i = 0; i < DIM; ++i)
        {
            box_[i]        = static_cast<real>(shape_[i]) * spacing_[i];
            invSpacing_[i] = 1.0 / spacing_[i];
        }

        const auto num_elements = shape_[XX] * shape_[YY] * shape_[ZZ];
        values_.resize(num_elements);
    }

private:
    //! Get the 1D index in the `values` array for a 3D grid position
    //!
    //! Throws `std::out_of_range` if 3D position is not within grid.
    size_t gridPositionToIndex(const size_t ix, const size_t iy, const size_t iz) const
    {
        const auto nx = static_cast<size_t>(shape_[XX]);
        const auto ny = static_cast<size_t>(shape_[YY]);
        const auto nz = static_cast<size_t>(shape_[ZZ]);

        if ((ix >= nx) || (iy >= ny) || (iz >= nz))
        {
            char buf[STRLEN];
            snprintf(buf,
                     STRLEN,
                     "Grid3d::_index: position (%lu, %lu, %lu) not within "
                     "grid of size (%lu, %lu, %lu)",
                     ix,
                     iy,
                     iz,
                     nx,
                     ny,
                     nz);

            throw std::out_of_range(buf);
        }

        return iz + (iy * nz) + (ix * (ny * nz));
    }

    //! Number of cells along each axis
    IVec shape_;

    //! Size of cells along each axis
    RVec spacing_;

    //! Position of grid corner in system
    RVec origin_ = { 0.0, 0.0, 0.0 };

    //! Physical size of grid
    RVec box_;

    //! 1.0 / spacing along each axis, used to compute indexing from position
    RVec invSpacing_;

    //! Values inside the 3D grid
    //!
    //! Stored in Z-Y-X order.
    std::vector<T> values_;
};

} // namespace flow
} // namespace gmx

#endif // MD_FLOW_FIELD_GRID
