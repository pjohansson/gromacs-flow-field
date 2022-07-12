#include <cstdio>

#include "grid.h"

template<typename T>
const T& flow::Grid3d<T>::at(const int ix,
                             const int iy,
                             const int iz) const
{
    return values.at(_index(
        static_cast<size_t>(ix),
        static_cast<size_t>(iy),
        static_cast<size_t>(iz)
    ));
}

template<typename T>
T& flow::Grid3d<T>::at(const int ix,
                       const int iy,
                       const int iz)
{
    return values.at(_index(
        static_cast<size_t>(ix),
        static_cast<size_t>(iy),
        static_cast<size_t>(iz)
    ));
}

static int
get_pos_in_grid_saturated(const size_t dim,
                          const rvec   r,
                          const ivec   shape,
                          const rvec   origin,
                          const rvec   spacing)
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
    const auto ix = get_pos_in_grid_saturated(XX, r, shape, origin, spacing);
    const auto iy = get_pos_in_grid_saturated(YY, r, shape, origin, spacing);
    const auto iz = get_pos_in_grid_saturated(ZZ, r, shape, origin, spacing);

    return at(ix, iy, iz);
}

template<typename T>
T& flow::Grid3d<T>::at_pos(const rvec r)
{
    const auto ix = get_pos_in_grid_saturated(XX, r, shape, origin, spacing);
    const auto iy = get_pos_in_grid_saturated(YY, r, shape, origin, spacing);
    const auto iz = get_pos_in_grid_saturated(ZZ, r, shape, origin, spacing);

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
template struct flow::Grid3d<float>;
template struct flow::Grid3d<double>;
