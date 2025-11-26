#include "grid.h"

#include <cmath>
#include <cstdio>

#include "flow_field.h"

namespace gmx
{
namespace flow
{

template<typename T>
const T& flow::Grid3d<T>::at(const int ix, const int iy, const int iz) const
{
    return values_.at(gridPositionToIndex(
            static_cast<size_t>(ix), static_cast<size_t>(iy), static_cast<size_t>(iz)));
}

template<typename T>
T& flow::Grid3d<T>::at(const int ix, const int iy, const int iz)
{
    return values_.at(gridPositionToIndex(
            static_cast<size_t>(ix), static_cast<size_t>(iy), static_cast<size_t>(iz)));
}

static int get_pos_in_grid_saturated(const size_t dim, const rvec r, const ivec shape, const rvec origin, const rvec spacing)
{
    auto i = static_cast<int>((r[dim] - origin[dim]) / spacing[dim]);

    if (i < 0)
    {
        i = 0;
    }

    if (i >= shape[dim])
    {
        i = shape[dim] - 1;
    }

    return i;
}

template<typename T>
const T& flow::Grid3d<T>::at_pos(const rvec r) const
{
    const auto ix = get_pos_in_grid_saturated(XX, r, shape_, origin_, spacing_);
    const auto iy = get_pos_in_grid_saturated(YY, r, shape_, origin_, spacing_);
    const auto iz = get_pos_in_grid_saturated(ZZ, r, shape_, origin_, spacing_);

    return at(ix, iy, iz);
}

template<typename T>
T& flow::Grid3d<T>::at_pos(const rvec r)
{
    const auto ix = get_pos_in_grid_saturated(XX, r, shape_, origin_, spacing_);
    const auto iy = get_pos_in_grid_saturated(YY, r, shape_, origin_, spacing_);
    const auto iz = get_pos_in_grid_saturated(ZZ, r, shape_, origin_, spacing_);

    return at(ix, iy, iz);
}

template<typename T>
T& flow::Grid3d<T>::at_pos_pbc(const rvec r0, const matrix box)
{
    rvec r_pbc = { std::fmod(r0[XX], box[XX][XX]),
                   std::fmod(r0[YY], box[YY][YY]),
                   std::fmod(r0[ZZ], box[ZZ][ZZ]) };

    for (size_t i = 0; i < DIM; ++i)
    {
        while (r_pbc[i] < 0.0)
        {
            r_pbc[i] += box[i][i];
        }
    }

    const auto ix = get_pos_in_grid_saturated(XX, r_pbc, shape_, origin_, spacing_);
    const auto iy = get_pos_in_grid_saturated(YY, r_pbc, shape_, origin_, spacing_);
    const auto iz = get_pos_in_grid_saturated(ZZ, r_pbc, shape_, origin_, spacing_);

    return at(ix, iy, iz);
}

template<typename T>
const T& flow::Grid3d<T>::at_pos_pbc(const rvec r0, const matrix box) const
{
    rvec r_pbc = { std::fmod(r0[XX], box[XX][XX]),
                   std::fmod(r0[YY], box[YY][YY]),
                   std::fmod(r0[ZZ], box[ZZ][ZZ]) };

    for (size_t i = 0; i < DIM; ++i)
    {
        while (r_pbc[i] < 0.0)
        {
            r_pbc[i] += box[i][i];
        }
    }

    const auto ix = get_pos_in_grid_saturated(XX, r_pbc, shape_, origin_, spacing_);
    const auto iy = get_pos_in_grid_saturated(YY, r_pbc, shape_, origin_, spacing_);
    const auto iz = get_pos_in_grid_saturated(ZZ, r_pbc, shape_, origin_, spacing_);

    return at(ix, iy, iz);
}


// Since this file does not know which types `T` to generate
// code for, we need to instantiate the structure for all
// needed types during compilation.
//
// Pros: All code is generated once and linked to. If this
//       code was inlined into the class definition (in grid.h)
//       each created object would have a copy of the code,
//       increasing their size.
//
// Cons: We need to know all the types in advance and instantiate
//       them as below.
//
// Thoughts: After finishing development, consider moving to
//           inlined methods?
//
// Also, this would not be an issue in Rust. God bless
// the borrow checker.
template class Grid3d<float>;
template class Grid3d<double>;
template class Grid3d<Bin>;

} // namespace flow
} // namespace gmx
