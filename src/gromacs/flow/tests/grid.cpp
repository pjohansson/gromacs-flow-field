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

/*! \brief Tests for Grid3d class routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gromacs/flow/grid.h"

#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{
namespace test
{
namespace
{

using namespace flow;

TEST(FlowGridTest, ShapeAndSpacingAresSet)
{
    const IVec   shape   = { 1, 2, 3 };
    const RVec   spacing = { 4.0, 5.0, 6.0 };
    const Grid3d grid(shape, spacing);

    EXPECT_EQ(grid.shape(), shape);
    EXPECT_EQ(grid.spacing(), spacing);
}

TEST(FlowGridTest, SetsCorrectBoxSize)
{
    const IVec   shape   = { 1, 2, 3 };
    const RVec   spacing = { 4.0, 5.0, 6.0 };
    const Grid3d grid(shape, spacing);

    EXPECT_FLOAT_EQ(grid.box()[XX], shape[XX] * spacing[XX]);
    EXPECT_FLOAT_EQ(grid.box()[YY], shape[YY] * spacing[YY]);
    EXPECT_FLOAT_EQ(grid.box()[ZZ], shape[ZZ] * spacing[ZZ]);
}

TEST(FlowGridTest, InitializesValuesToCorrectNumValues)
{
    constexpr int nx = 3, ny = 5, nz = 7;
    const Grid3d  grid({ nx, ny, nz }, { -1.0, -1.0, -1.0 });
    EXPECT_EQ(nx * ny * nz, grid.values().size());
}

TEST(FlowGridTest, GridShapeMustBePositiveInAllDirections)
{
    EXPECT_THROW(Grid3d({ 0, 1, 1 }, { 1.0, 1.0, 1.0 }), std::invalid_argument);
    EXPECT_THROW(Grid3d({ 1, 0, 1 }, { 1.0, 1.0, 1.0 }), std::invalid_argument);
    EXPECT_THROW(Grid3d({ 1, 1, 0 }, { 1.0, 1.0, 1.0 }), std::invalid_argument);

    EXPECT_THROW(Grid3d({ -1, 1, 1 }, { 1.0, 1.0, 1.0 }), std::invalid_argument);
    EXPECT_THROW(Grid3d({ 1, -1, 1 }, { 1.0, 1.0, 1.0 }), std::invalid_argument);
    EXPECT_THROW(Grid3d({ 1, 1, -1 }, { 1.0, 1.0, 1.0 }), std::invalid_argument);
}

TEST(FlowGridTest, BinVolumeCalc)
{
    const IVec   shape   = { 3, 5, 7 };
    const RVec   spacing = { 11.0, 13.0, 17.0 };
    const Grid3d grid(shape, spacing);

    const double binVolume = spacing[XX] * spacing[YY] * spacing[ZZ];
    EXPECT_FLOAT_EQ(binVolume, grid.binVolume());
}

TEST(FlowGridTest, AssignSetsAllCellsToValue)
{
    constexpr double value = 5.0;

    Grid3d grid({ 3, 5, 7 }, { -1.0, -1.0, -1.0 });
    grid.assign(value);

    for (const double v : grid.values())
    {
        EXPECT_FLOAT_EQ(v, value);
    }
}

TEST(FlowGridTest, DataInArrayUsesOrderingZYX)
{
    const size_t nx = 3, ny = 5, nz = 7;

    Grid3d grid({ nx, ny, nz }, { -1.0, -1.0, -1.0 });

    std::iota(grid.values().begin(), grid.values().end(), 0.0);

    // Z increases first
    EXPECT_FLOAT_EQ(0.0, grid.at(0, 0, 0));
    EXPECT_FLOAT_EQ(1.0, grid.at(0, 0, 1));
    EXPECT_FLOAT_EQ(2.0, grid.at(0, 0, 2));

    // ... then Y (starting at =nz)
    EXPECT_FLOAT_EQ(static_cast<double>(nz), grid.at(0, 1, 0));

    // ... finally x (starting at =ny*nz)
    EXPECT_FLOAT_EQ(static_cast<double>(ny * nz), grid.at(1, 0, 0));
}

TEST(FlowGridTest, OutOfBoundsAccessThrows)
{
    const size_t nx = 3, ny = 5, nz = 7;

    Grid3d grid({ nx, ny, nz }, { -1.0, -1.0, -1.0 });

    EXPECT_THROW(grid.at(nx, 0, 0), std::out_of_range);
    EXPECT_THROW(grid.at(0, ny, 0), std::out_of_range);
    EXPECT_THROW(grid.at(0, 0, nz), std::out_of_range);
}

TEST(FlowGridTest, ContainsWorks)
{
    const IVec shape{ 3, 5, 7 };
    const RVec box{ 5.0, 7.0, 11.0 };
    const RVec spacing{ box[XX] / static_cast<real>(shape[XX]),
                        box[YY] / static_cast<real>(shape[YY]),
                        box[ZZ] / static_cast<real>(shape[ZZ]) };

    const Grid3d grid(shape, spacing);

    EXPECT_TRUE(grid.contains(RVec{ 0.0f, 0.0f, 0.0f }));
    EXPECT_TRUE(grid.contains(box));

    constexpr real eps = 1e-3; // small adjustment to move just out of box
    EXPECT_FALSE(grid.contains(RVec{ -eps, 0.0, 0.0 }));
    EXPECT_FALSE(grid.contains(RVec{ 0.0, -eps, 0.0 }));
    EXPECT_FALSE(grid.contains(RVec{ 0.0, 0.0, -eps }));
    EXPECT_FALSE(grid.contains(RVec{ box[XX] + eps, 0.0, 0.0 }));
    EXPECT_FALSE(grid.contains(RVec{ 0.0, box[YY] + eps, 0.0 }));
    EXPECT_FALSE(grid.contains(RVec{ 0.0, 0.0, box[ZZ] + eps }));
}

TEST(FlowGridTest, AtPositionWorks)
{
    const IVec shape   = { 7, 11, 13 };
    const RVec spacing = { 0.5, 1.0, 2.0 };
    Grid3d     grid(shape, spacing);
    std::iota(grid.values().begin(), grid.values().end(), 0.0);

    // Grid positions to test: must lie within shape defined above
    // (outside values checked in other test below)
    const IVec gridPositionsToTest[] = {
        // Edge cases: at corners
        { 0, 0, 0 },
        { 6, 0, 0 },
        { 0, 10, 0 },
        { 0, 0, 12 },
        { 6, 10, 12 },
        // Some internal positions
        { 1, 1, 1 },
        { 3, 5, 6 },
        { 5, 0, 11 },
        { 0, 9, 11 },
        { 5, 9, 0 },
        { 5, 9, 11 },
    };

    for (const IVec& gridPosition : gridPositionsToTest)
    {
        // Get position in center of chosen bin
        const RVec position = { (static_cast<real>(gridPosition[XX]) + 0.5f) * spacing[XX],
                                (static_cast<real>(gridPosition[YY]) + 0.5f) * spacing[YY],
                                (static_cast<real>(gridPosition[ZZ]) + 0.5f) * spacing[ZZ] };

        const int index = gridPosition[ZZ] + (gridPosition[YY] * shape[ZZ])
                          + (gridPosition[XX] * shape[YY] * shape[ZZ]);
        EXPECT_EQ(grid.atPosition(position), grid.values()[index]);
    }
}

TEST(FlowGridTest, AtPositionPBCPutsPosInBox)
{
    const IVec shape   = { 7, 11, 13 };
    const RVec spacing = { 0.5, 1.0, 2.0 };
    Grid3d     grid(shape, spacing);
    std::iota(grid.values().begin(), grid.values().end(), 0.0);

    const RVec& box = grid.box();
    const real  x0  = 1.0 * spacing[XX];
    const real  y0  = 2.0 * spacing[YY];
    const real  z0  = 3.0 * spacing[ZZ];

    // Along x
    EXPECT_FLOAT_EQ(grid.atPosition(RVec{ x0, y0, z0 }), grid.atPosition(RVec{ x0 + box[XX], y0, z0 }));
    EXPECT_FLOAT_EQ(grid.atPosition(RVec{ x0, y0, z0 }), grid.atPosition(RVec{ x0 - box[XX], y0, z0 }));

    // Along y
    EXPECT_FLOAT_EQ(grid.atPosition(RVec{ x0, y0, z0 }), grid.atPosition(RVec{ x0, y0 + box[YY], z0 }));
    EXPECT_FLOAT_EQ(grid.atPosition(RVec{ x0, y0, z0 }), grid.atPosition(RVec{ x0, y0 - box[YY], z0 }));

    // Along z
    EXPECT_FLOAT_EQ(grid.atPosition(RVec{ x0, y0, z0 }), grid.atPosition(RVec{ x0, y0, z0 + box[ZZ] }));
    EXPECT_FLOAT_EQ(grid.atPosition(RVec{ x0, y0, z0 }), grid.atPosition(RVec{ x0, y0, z0 - box[ZZ] }));

    // Along all dimensions
    EXPECT_FLOAT_EQ(grid.atPosition(RVec{ x0, y0, z0 }),
                    grid.atPosition(RVec{ x0 + box[XX], y0 + box[YY], z0 + box[ZZ] }));
    EXPECT_FLOAT_EQ(grid.atPosition(RVec{ x0, y0, z0 }),
                    grid.atPosition(RVec{ x0 - box[XX], y0 - box[YY], z0 - box[ZZ] }));

    // // Along all dimensions, multiple shifts
    EXPECT_FLOAT_EQ(
            grid.atPosition(RVec{ x0, y0, z0 }),
            grid.atPosition(RVec{ x0 + (3.0f * box[XX]), y0 + (5.0f * box[YY]), z0 + (7.0f * box[ZZ]) }));
    EXPECT_FLOAT_EQ(
            grid.atPosition(RVec{ x0, y0, z0 }),
            grid.atPosition(RVec{ x0 - (7.0f * box[XX]), y0 - (3.0f * box[YY]), z0 - (5.0f * box[ZZ]) }));
}

TEST(FlowGridTest, IndexFromPositionWorks)
{
    const IVec shape   = { 7, 11, 13 };
    const RVec spacing = { 0.5, 1.0, 2.0 };
    Grid3d     grid(shape, spacing);
    std::iota(grid.values().begin(), grid.values().end(), 0.0);

    // Grid positions to test: must lie within shape defined above
    // (outside values checked in other test below)
    const IVec gridPositionsToTest[] = {
        // Edge cases: at corners
        { 0, 0, 0 },
        { 6, 0, 0 },
        { 0, 10, 0 },
        { 0, 0, 12 },
        { 6, 10, 12 },
        // Some internal positions
        { 1, 1, 1 },
        { 3, 5, 6 },
        { 5, 0, 11 },
        { 0, 9, 11 },
        { 5, 9, 0 },
        { 5, 9, 11 },
        // Some outside positions (PBC adjusted to put inside box)
        { -1, -1, -1 },
        { 57, 103, 85 },
        { -57, 103, 85 },
    };

    for (const IVec& gridPosition : gridPositionsToTest)
    {
        // Get position in center of chosen bin
        const RVec position = { (static_cast<real>(gridPosition[XX]) + 0.5f) * spacing[XX],
                                (static_cast<real>(gridPosition[YY]) + 0.5f) * spacing[YY],
                                (static_cast<real>(gridPosition[ZZ]) + 0.5f) * spacing[ZZ] };

        int ix = gridPosition[XX] % shape[XX];
        int iy = gridPosition[YY] % shape[YY];
        int iz = gridPosition[ZZ] % shape[ZZ];

        while (ix < 0)
        {
            ix += shape[XX];
        }
        while (iy < 0)
        {
            iy += shape[YY];
        }
        while (iz < 0)
        {
            iz += shape[ZZ];
        }

        const int index = iz + (iy * shape[ZZ]) + (ix * shape[YY] * shape[ZZ]);
        EXPECT_EQ(grid.indexFromPosition(position), index);
    }
}

TEST(FlowGridTest, SetBoxUpdatesSpacings)
{
    const IVec shape   = { 7, 11, 13 };
    const RVec spacing = { 0.5, 1.0, 2.0 };
    Grid3d     grid(shape, spacing);

    const matrix newBox = { { 2.0f * shape[XX] * spacing[XX], 0.0, 0.0 },
                            { 0.0, 2.0f * shape[YY] * spacing[YY], 0.0 },
                            { 0.0, 0.0, 2.0f * shape[ZZ] * spacing[ZZ] } };
    grid.setBox(newBox);

    for (int i = 0; i < DIM; ++i)
    {
        EXPECT_EQ(grid.box()[i], newBox[i][i]);
    }
    EXPECT_EQ(grid.shape(), shape);
    EXPECT_EQ(grid.spacing(), 2.0f * spacing);
}

} // namespace
} // namespace test
} // namespace gmx
