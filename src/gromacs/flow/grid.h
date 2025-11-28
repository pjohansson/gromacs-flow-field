#ifndef GMX_FLOW_GRID_H
#define GMX_FLOW_GRID_H

#include <vector>

#include "gromacs/math/vectypes.h"

namespace gmx
{

template<typename T>
class ArrayRef;

namespace flow
{

//! \brief Class which manages values inside a 3D grid.
template<typename T = double>
class Grid3d
{
public:
    //! Empty constructor
    Grid3d() {}

    //! Default constructor of grid, does not set values
    Grid3d(const IVec& shape, const RVec& spacing);

    //! Assign \p value to all cells in the grid.
    void assign(const T& value);

    //! Grid cell reference accessors.
    //!
    //! Throws `std::out_of_range` if 3D position is not within grid.
    const T& at(int ix, int iy, int iz) const;
    T&       at(int ix, int iy, int iz);

    //! Grid cell reference accessors from particle \p position.
    // Applies PBC correction to put \p inside the \c box_.
    const T& atPosition(const RVec& position) const;
    T&       atPosition(const RVec& position);

    //! Return the 1d index of \c values() corresponding to input \p position.
    size_t indexFromPosition(const RVec& position) const;

    //! Check whether a position is contained within the grid
    bool contains(const RVec& position) const;

    //! Return the bin volume
    real binVolume() const noexcept;

    //! Return the grid shape.
    const IVec& shape() const { return shape_; }

    //! Return the grid spacing.
    const RVec& spacing() const { return spacing_; }

    //! Return the box containing the grid.
    const RVec& box() const { return box_; }

    //! Return a reference to the stored values.
    ArrayRef<const T> values() const;
    ArrayRef<T>       values();

private:
    //! Get the 1D index in the `values` array for a 3D grid position
    //!
    //! Throws `std::out_of_range` if 3D position is not within grid.
    size_t gridPositionToIndex(size_t ix, size_t iy, size_t iz) const;

    //! Number of cells along each axis
    IVec shape_;

    //! Size of cells along each axis
    RVec spacing_;

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

#endif // GMX_FLOW_GRID_H
