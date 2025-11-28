#include "grid.h"

#include <cmath>
#include <cstdio>

#include <stdexcept>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/stringutil.h"

#include "flow_field.h"

namespace gmx
{
namespace flow
{

template<typename T>
Grid3d<T>::Grid3d(const IVec& shape, const RVec& spacing) : shape_{ shape }, spacing_{ spacing }
{
    if ((shape_[XX] < 1) || (shape_[YY] < 1) || (shape_[ZZ] < 1))
    {
        throw std::invalid_argument(formatString(
                "Grid3d::Grid3d: shape (%d, %d, %d) must be positive along all directions",
                shape_[XX],
                shape_[YY],
                shape_[ZZ]));
    }

    for (int i = 0; i < DIM; ++i)
    {
        box_[i]        = static_cast<real>(shape_[i]) * spacing_[i];
        invSpacing_[i] = 1.0 / spacing_[i];
    }

    values_.resize(shape_[XX] * shape_[YY] * shape_[ZZ]);
}

template<typename T>
void Grid3d<T>::assign(const T& value)
{
    for (auto& v : values_)
    {
        v = value;
    }
}

template<typename T>
const T& Grid3d<T>::at(const int ix, const int iy, const int iz) const
{
    return values_.at(gridPositionToIndex(
            static_cast<size_t>(ix), static_cast<size_t>(iy), static_cast<size_t>(iz)));
}

template<typename T>
T& Grid3d<T>::at(const int ix, const int iy, const int iz)
{
    return const_cast<T&>(const_cast<const Grid3d*>(this)->at(ix, iy, iz));
}

template<typename T>
const T& Grid3d<T>::atPosition(const RVec& position) const
{

    return values_[indexFromPosition(position)];
}

template<typename T>
T& Grid3d<T>::atPosition(const RVec& position)
{
    return values_[indexFromPosition(position)];
}

template<typename T>
size_t Grid3d<T>::indexFromPosition(const RVec& position) const
{
    IVec gridPosition = { static_cast<int>(std::floor(position[XX] * invSpacing_[XX])) % shape_[XX],
                          static_cast<int>(std::floor(position[YY] * invSpacing_[YY])) % shape_[YY],
                          static_cast<int>(std::floor(position[ZZ] * invSpacing_[ZZ])) % shape_[ZZ] };

    for (int i = 0; i < DIM; ++i)
    {
        while (gridPosition[i] < 0)
        {
            gridPosition[i] += shape_[i];
        }
    }

    return gridPosition[ZZ] + (gridPosition[YY] * shape_[ZZ])
           + (gridPosition[XX] * shape_[YY] * shape_[ZZ]);
}
template<typename T>
real Grid3d<T>::binVolume() const noexcept
{
    return spacing_[XX] * spacing_[YY] * spacing_[ZZ];
}

template<typename T>
bool Grid3d<T>::contains(const RVec& r) const
{
    return (r[XX] >= 0.0 && r[XX] <= box_[XX] && r[YY] >= 0.0 && r[YY] <= box_[YY] && r[ZZ] >= 0.0
            && r[ZZ] <= box_[ZZ]);
}

template<typename T>
ArrayRef<T> Grid3d<T>::values()
{
    return values_;
}

template<typename T>
ArrayRef<const T> Grid3d<T>::values() const
{
    return values_;
}

template<typename T>
size_t Grid3d<T>::gridPositionToIndex(const size_t ix, const size_t iy, const size_t iz) const
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

    return iz + (iy * nz) + (ix * ny * nz);
}

template class Grid3d<float>;
template class Grid3d<double>;
template class Grid3d<Bin>;

} // namespace flow
} // namespace gmx
