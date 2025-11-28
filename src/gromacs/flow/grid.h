/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \brief Declarations of Grid3d class routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

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
