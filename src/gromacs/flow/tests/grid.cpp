#include "gromacs/flow/grid.h"

#include <optional>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

namespace gmx
{
namespace test
{
namespace
{

using namespace gmx;
using namespace flow;

TEST(FlowGridTest, InitializesWithInput)
{
    const auto grid = Grid3d({ 1, 2, 3 }, { 4.0, 5.0, 6.0 }, {});

    EXPECT_EQ(1, grid.shape()[XX]);
    EXPECT_EQ(2, grid.shape()[YY]);
    EXPECT_EQ(3, grid.shape()[ZZ]);

    EXPECT_FLOAT_EQ(4.0, grid.spacing()[XX]);
    EXPECT_FLOAT_EQ(5.0, grid.spacing()[YY]);
    EXPECT_FLOAT_EQ(6.0, grid.spacing()[ZZ]);
}

TEST(FlowGridTest, InitializesWithInverseBinSpacing)
{
    const auto grid = Grid3d({ 1, 2, 3 }, { 4.0, 5.0, 6.0 }, {});

    EXPECT_FLOAT_EQ(1.0 / 4.0, grid.invSpacing()[XX]);
    EXPECT_FLOAT_EQ(1.0 / 5.0, grid.invSpacing()[YY]);
    EXPECT_FLOAT_EQ(1.0 / 6.0, grid.invSpacing()[ZZ]);
}

TEST(FlowGridTest, InitializesWithCorrectGridSize)
{
    const int nx = 3, ny = 5, nz = 7;

    const auto grid = Grid3d({ nx, ny, nz }, { 4.0, 5.0, 6.0 }, {});
    EXPECT_EQ(nx * ny * nz, grid.values().size());
}

TEST(FlowGridTest, GridShapeMustBePositiveInAllDirections)
{
    const RVec spacing = { 1.0, 1.0, 1.0 };

    EXPECT_THROW(Grid3d({ 0, 1, 1 }, spacing, {}), std::invalid_argument);
    EXPECT_THROW(Grid3d({ 1, 0, 1 }, spacing, {}), std::invalid_argument);
    EXPECT_THROW(Grid3d({ 1, 1, 0 }, spacing, {}), std::invalid_argument);

    EXPECT_THROW(Grid3d({ -1, 1, 1 }, spacing, {}), std::invalid_argument);
    EXPECT_THROW(Grid3d({ 1, -1, 1 }, spacing, {}), std::invalid_argument);
    EXPECT_THROW(Grid3d({ 1, 1, -1 }, spacing, {}), std::invalid_argument);
}

TEST(FlowGridTest, GridOriginIsZeroByDefault)
{
    const IVec shape   = { 3, 5, 7 };
    const RVec spacing = { 11.0, 13.0, 17.0 };

    const auto grid_default = Grid3d(shape, spacing, {});
    EXPECT_FLOAT_EQ(grid_default.origin()[XX], 0.0);
    EXPECT_FLOAT_EQ(grid_default.origin()[YY], 0.0);
    EXPECT_FLOAT_EQ(grid_default.origin()[ZZ], 0.0);
}

TEST(FlowGridTest, GridOriginFromOptionalArgument)
{
    const IVec shape   = { 3, 5, 7 };
    const RVec spacing = { 11.0, 13.0, 17.0 };

    const real x0 = 19.0, y0 = 23.0, z0 = 29.0;

    const auto grid = Grid3d(shape, spacing, RVec{ x0, y0, z0 });
    EXPECT_FLOAT_EQ(grid.origin()[XX], x0);
    EXPECT_FLOAT_EQ(grid.origin()[YY], y0);
    EXPECT_FLOAT_EQ(grid.origin()[ZZ], z0);
}

TEST(FlowGridTest, BinVolumeCalc)
{
    const IVec shape   = { 3, 5, 7 };
    const RVec spacing = { 11.0, 13.0, 17.0 };

    const auto grid = Grid3d(shape, spacing, {});

    const auto bin_volume = spacing[XX] * spacing[YY] * spacing[ZZ];
    EXPECT_FLOAT_EQ(bin_volume, grid.bin_volume());
}

TEST(FlowGridTest, AssignSetsAllCellsToValue)
{
    const size_t nx = 3, ny = 5, nz = 7;

    auto grid = Grid3d({ nx, ny, nz }, { 4.0, 5.0, 6.0 }, {});

    const double value = 5.0;
    grid.assign(value);

    for (const auto& v : grid.values())
    {
        EXPECT_FLOAT_EQ(v, value);
    }
}

TEST(FlowGridTest, DataInArrayUsesOrderingZYX)
{
    const size_t nx = 3, ny = 5, nz = 7;

    auto grid = Grid3d({ nx, ny, nz }, { 4.0, 5.0, 6.0 }, {});

    double val = 0.0;
    for (auto& v : grid.values())
    {
        v = val;
        val += 1.0;
    }

    // Z increases first
    EXPECT_FLOAT_EQ(0.0, grid.at(0, 0, 0));
    EXPECT_FLOAT_EQ(1.0, grid.at(0, 0, 1));
    EXPECT_FLOAT_EQ(2.0, grid.at(0, 0, 2));

    // ... then Y (starting at =nz)
    EXPECT_FLOAT_EQ(static_cast<double>(nz), grid.at(0, 1, 0));

    // ... finally x (starting at =ny*nz)
    EXPECT_FLOAT_EQ(static_cast<double>(ny * nz), grid.at(1, 0, 0));
}

TEST(FlowGridTest, AtMethodCanGetValue)
{
    const size_t nx = 3, ny = 5, nz = 7;

    auto grid = Grid3d({ nx, ny, nz }, { 4.0, 5.0, 6.0 }, {});
    grid.assign(0.0);

    grid.values().front() = 3.0;
    EXPECT_FLOAT_EQ(3.0, grid.at(0, 0, 0));

    grid.values().back() = 5.0;
    EXPECT_FLOAT_EQ(5.0, grid.at(nx - 1, ny - 1, nz - 1));
}

TEST(FlowGridTest, AtMethodCanSetValue)
{
    const size_t nx = 3, ny = 5, nz = 7;

    auto grid = Grid3d({ nx, ny, nz }, { 4.0, 5.0, 6.0 }, {});
    grid.assign(0.0);

    const double value = 15.0;
    const size_t ix = 1, iy = 2, iz = 3;

    grid.at(ix, iy, iz) = value;
    EXPECT_FLOAT_EQ(value, grid.at(ix, iy, iz));
}

TEST(FlowGridTest, AtMethodCanGetReference)
{
    const size_t nx = 3, ny = 5, nz = 7;

    auto grid = Grid3d({ nx, ny, nz }, { 4.0, 5.0, 6.0 }, {});
    grid.assign(0.0);

    const double value = 15.0;
    const size_t ix = 1, iy = 2, iz = 3;

    const double& ref   = grid.at(ix, iy, iz);
    grid.at(ix, iy, iz) = value;

    EXPECT_FLOAT_EQ(value, static_cast<double>(ref));
}

TEST(FlowGridTest, AtMethodCanSetValueToReference)
{
    const size_t nx = 3, ny = 5, nz = 7;

    auto grid = Grid3d({ nx, ny, nz }, { 4.0, 5.0, 6.0 }, {});

    const double value = 15.0;
    const size_t ix = 1, iy = 2, iz = 3;

    auto& ref = grid.at(ix, iy, iz);
    ref       = value;

    EXPECT_FLOAT_EQ(value, grid.at(ix, iy, iz));
}

TEST(FlowGridTest, AtMethodWithAutoReturnsValueNotReference)
{
    const size_t nx = 3, ny = 5, nz = 7;

    auto grid = Grid3d({ nx, ny, nz }, { 4.0, 5.0, 6.0 }, {});

    const size_t ix = 1, iy = 2, iz = 3;

    const double value1 = 15.0;
    const double value2 = 2.0 * value1;
    const double value3 = 2.0 * value2;

    grid.at(ix, iy, iz) = value1;
    auto non_ref        = grid.at(ix, iy, iz);

    grid.at(ix, iy, iz) = value2;
    auto& ref           = grid.at(ix, iy, iz);

    grid.at(ix, iy, iz) = value3;

    EXPECT_FLOAT_EQ(value1, non_ref);
    EXPECT_FLOAT_EQ(value3, ref);
}

TEST(FlowGridTest, OutOfBoundsAccessThrows)
{
    const size_t nx = 3, ny = 5, nz = 7;

    auto grid = Grid3d({ nx, ny, nz }, { 4.0, 5.0, 6.0 }, {});

    EXPECT_THROW(grid.at(nx, 0, 0), std::out_of_range);
    EXPECT_THROW(grid.at(0, ny, 0), std::out_of_range);
    EXPECT_THROW(grid.at(0, 0, nz), std::out_of_range);
}

TEST(FlowGridTest, ContainsWithRawPointer)
{
    const size_t nx = 3, ny = 5, nz = 7;

    const real box_x = 5.0, box_y = 7.0, box_z = 11.0;

    const real x0 = 1.0, y0 = 2.0, z0 = 3.0;
    const RVec origin{ x0, y0, z0 };

    const real dx = box_x / static_cast<real>(nx), dy = box_y / static_cast<real>(ny),
               dz = box_z / static_cast<real>(nz);

    const auto grid = Grid3d({ nx, ny, nz }, { dx, dy, dz }, origin);

    constexpr real d     = 1e-3;
    const rvec     r_in1 = { x0 + d, y0 + d, z0 + d },
               r_in2     = { x0 + box_x - d, y0 + box_y - d, z0 + box_z - d },
               r_out1 = { x0 - d, y0 + d, z0 + d }, r_out2 = { x0 + d, y0 - d, z0 + d },
               r_out3 = { x0 + d, y0 + d, z0 - d },
               r_out4 = { x0 + box_x + d, y0 + box_y - d, z0 + box_z - d },
               r_out5 = { x0 + box_x - d, y0 + box_y + d, z0 + box_z - d },
               r_out6 = { x0 + box_x - d, y0 + box_y - d, z0 + box_z + d };

    EXPECT_TRUE(grid.contains(r_in1));
    EXPECT_TRUE(grid.contains(r_in2));

    EXPECT_FALSE(grid.contains(r_out1));
    EXPECT_FALSE(grid.contains(r_out2));
    EXPECT_FALSE(grid.contains(r_out3));
    EXPECT_FALSE(grid.contains(r_out4));
    EXPECT_FALSE(grid.contains(r_out5));
    EXPECT_FALSE(grid.contains(r_out6));
}

TEST(FlowGridTest, ContainsWithRVec)
{
    const size_t nx = 3, ny = 5, nz = 7;

    const real box_x = 5.0, box_y = 7.0, box_z = 11.0;

    const real x0 = 1.0, y0 = 2.0, z0 = 3.0;
    const RVec origin{ x0, y0, z0 };

    const real dx = box_x / static_cast<real>(nx), dy = box_y / static_cast<real>(ny),
               dz = box_z / static_cast<real>(nz);

    const auto grid = Grid3d({ nx, ny, nz }, { dx, dy, dz }, origin);

    constexpr real d     = 1e-3;
    const RVec     r_in1 = { x0 + d, y0 + d, z0 + d },
               r_in2     = { x0 + box_x - d, y0 + box_y - d, z0 + box_z - d },
               r_out1 = { x0 - d, y0 + d, z0 + d }, r_out2 = { x0 + d, y0 - d, z0 + d },
               r_out3 = { x0 + d, y0 + d, z0 - d },
               r_out4 = { x0 + box_x + d, y0 + box_y - d, z0 + box_z - d },
               r_out5 = { x0 + box_x - d, y0 + box_y + d, z0 + box_z - d },
               r_out6 = { x0 + box_x - d, y0 + box_y - d, z0 + box_z + d };

    EXPECT_TRUE(grid.contains(r_in1));
    EXPECT_TRUE(grid.contains(r_in2));

    EXPECT_FALSE(grid.contains(r_out1));
    EXPECT_FALSE(grid.contains(r_out2));
    EXPECT_FALSE(grid.contains(r_out3));
    EXPECT_FALSE(grid.contains(r_out4));
    EXPECT_FALSE(grid.contains(r_out5));
    EXPECT_FALSE(grid.contains(r_out6));
}

TEST(FlowGridTest, AtPositionInsideGrid)
{
    const int nx = 4, ny = 4, nz = 4;

    const real dx = 0.5, dy = 1.0, dz = 2.0;

    auto grid = Grid3d({ nx, ny, nz }, { dx, dy, dz }, {});

    double val = 0.0;
    for (auto& v : grid.values())
    {
        v = val;
        val += 1.0;
    }

    EXPECT_FLOAT_EQ(grid.at(0, 0, 0), grid.at_pos(RVec{ 0.25, 0.5, 1.0 }));
    EXPECT_FLOAT_EQ(grid.at(1, 0, 0), grid.at_pos(RVec{ 0.75, 0.5, 1.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 0, 0), grid.at_pos(RVec{ 1.25, 0.5, 1.0 }));
    EXPECT_FLOAT_EQ(grid.at(3, 0, 0), grid.at_pos(RVec{ 1.75, 0.5, 1.0 }));

    EXPECT_FLOAT_EQ(grid.at(2, 1, 0), grid.at_pos(RVec{ 1.25, 1.5, 1.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 2, 0), grid.at_pos(RVec{ 1.25, 2.5, 1.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 3, 0), grid.at_pos(RVec{ 1.25, 3.5, 1.0 }));

    EXPECT_FLOAT_EQ(grid.at(2, 2, 1), grid.at_pos(RVec{ 1.25, 2.5, 3.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 2, 2), grid.at_pos(RVec{ 1.25, 2.5, 5.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 2, 3), grid.at_pos(RVec{ 1.25, 2.5, 7.0 }));
}

TEST(FlowGridTest, AtPositionOutsideGridSaturatesAtEdge)
{
    const int nx = 4, ny = 4, nz = 4;

    const real dx = 0.5, dy = 1.0, dz = 2.0;

    auto grid = Grid3d({ nx, ny, nz }, { dx, dy, dz }, {});

    double val = 0.0;
    for (auto& v : grid.values())
    {
        v = val;
        val += 1.0;
    }

    EXPECT_FLOAT_EQ(grid.at(0, 0, 0), grid.at_pos(RVec{ -0.25, 0.5, 1.0 }));
    EXPECT_FLOAT_EQ(grid.at(0, 0, 0), grid.at_pos(RVec{ 0.25, -0.5, 1.0 }));
    EXPECT_FLOAT_EQ(grid.at(0, 0, 0), grid.at_pos(RVec{ 0.25, 0.5, -1.0 }));

    EXPECT_FLOAT_EQ(grid.at(3, 0, 0), grid.at_pos(RVec{ 10.0, 0.5, 1.0 }));
    EXPECT_FLOAT_EQ(grid.at(0, 3, 0), grid.at_pos(RVec{ 0.25, 10.0, 1.0 }));
    EXPECT_FLOAT_EQ(grid.at(0, 0, 3), grid.at_pos(RVec{ 0.25, 0.5, 10.0 }));
}

TEST(FlowGridTest, AtPositionInsideGridWithShiftedOrigin)
{
    const int nx = 4, ny = 4, nz = 4;

    const real dx = 0.5, dy = 1.0, dz = 2.0;

    const real x0 = 10.0, y0 = 20.0, z0 = 30.0;

    auto grid = Grid3d({ nx, ny, nz }, { dx, dy, dz }, RVec{ x0, y0, z0 });

    double val = 0.0;
    for (auto& v : grid.values())
    {
        v = val;
        val += 1.0;
    }

    EXPECT_FLOAT_EQ(grid.at(0, 0, 0), grid.at_pos(RVec{ 10.25, 20.5, 31.0 }));
    EXPECT_FLOAT_EQ(grid.at(1, 0, 0), grid.at_pos(RVec{ 10.75, 20.5, 31.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 0, 0), grid.at_pos(RVec{ 11.25, 20.5, 31.0 }));
    EXPECT_FLOAT_EQ(grid.at(3, 0, 0), grid.at_pos(RVec{ 11.75, 20.5, 31.0 }));

    EXPECT_FLOAT_EQ(grid.at(2, 1, 0), grid.at_pos(RVec{ 11.25, 21.5, 31.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 2, 0), grid.at_pos(RVec{ 11.25, 22.5, 31.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 3, 0), grid.at_pos(RVec{ 11.25, 23.5, 31.0 }));

    EXPECT_FLOAT_EQ(grid.at(2, 2, 1), grid.at_pos(RVec{ 11.25, 22.5, 33.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 2, 2), grid.at_pos(RVec{ 11.25, 22.5, 35.0 }));
    EXPECT_FLOAT_EQ(grid.at(2, 2, 3), grid.at_pos(RVec{ 11.25, 22.5, 37.0 }));
}

TEST(FlowGridTest, AtPositionWorksWithConst)
{
    const int nx = 4, ny = 4, nz = 4;

    const real dx = 0.5, dy = 1.0, dz = 2.0;

    auto grid1 = Grid3d({ nx, ny, nz }, { dx, dy, dz }, {});

    double val = 0.0;
    for (auto& v : grid1.values())
    {
        v = val;
        val += 1.0;
    }

    const auto grid2 = grid1;

    EXPECT_FLOAT_EQ(grid2.at(3, 0, 0), grid2.at_pos(RVec{ 1.75, 0.5, 1.0 }));
}

TEST(FlowGridTest, AtPositionPBCPutsPosInBox)
{
    const int nx = 7, ny = 11, nz = 13;

    const real dx = 0.5, dy = 1.0, dz = 2.0;

    const real x0 = 10.0, y0 = 20.0, z0 = 30.0;

    const real box_x = 40.0, box_y = 50.0, box_z = 60.0;

    const matrix box = { { box_x, 0.0, 0.0 }, { 0.0, box_y, 0.0 }, { 0.0, 0.0, box_z } };

    auto grid = Grid3d({ nx, ny, nz }, { dx, dy, dz }, RVec{ x0, y0, z0 });

    double val = 0.0;
    for (auto& v : grid.values())
    {
        v = val;
        val += 1.0;
    }

    const real x = x0 + (1.0 * dx), y = y0 + (2.0 * dy), z = z0 + (1.0 * dz);

    // Along x
    EXPECT_FLOAT_EQ(grid.at_pos(RVec{ x, y, z }), grid.at_pos_pbc(RVec{ x + box_x, y, z }, box));
    EXPECT_FLOAT_EQ(grid.at_pos(RVec{ x, y, z }), grid.at_pos_pbc(RVec{ x - box_x, y, z }, box));

    // Along y
    EXPECT_FLOAT_EQ(grid.at_pos(RVec{ x, y, z }), grid.at_pos_pbc(RVec{ x, y + box_y, z }, box));
    EXPECT_FLOAT_EQ(grid.at_pos(RVec{ x, y, z }), grid.at_pos_pbc(RVec{ x, y - box_y, z }, box));

    // Along z
    EXPECT_FLOAT_EQ(grid.at_pos(RVec{ x, y, z }), grid.at_pos_pbc(RVec{ x, y, z + box_z }, box));
    EXPECT_FLOAT_EQ(grid.at_pos(RVec{ x, y, z }), grid.at_pos_pbc(RVec{ x, y, z - box_z }, box));

    // Along all dimensions
    EXPECT_FLOAT_EQ(grid.at_pos(RVec{ x, y, z }),
                    grid.at_pos_pbc(RVec{ x + box_x, y + box_y, z + box_z }, box));
    EXPECT_FLOAT_EQ(grid.at_pos(RVec{ x, y, z }),
                    grid.at_pos_pbc(RVec{ x - box_x, y - box_y, z - box_z }, box));

    // Along all dimensions, multiple shifts
    EXPECT_FLOAT_EQ(
            grid.at_pos(RVec{ x, y, z }),
            grid.at_pos_pbc(RVec{ x + (3.0f * box_x), y + (5.0f * box_y), z + (7.0f * box_z) }, box));
    EXPECT_FLOAT_EQ(
            grid.at_pos(RVec{ x, y, z }),
            grid.at_pos_pbc(RVec{ x - (7.0f * box_x), y - (3.0f * box_y), z - (5.0f * box_z) }, box));
}

TEST(FlowGridTest, AtPositionWithPBCWorksWithConst)
{
    const int nx = 4, ny = 4, nz = 4;

    const real dx = 0.5, dy = 1.0, dz = 2.0;

    const real box_x = 40.0, box_y = 50.0, box_z = 60.0;

    const matrix box = { { box_x, 0.0, 0.0 }, { 0.0, box_y, 0.0 }, { 0.0, 0.0, box_z } };

    auto grid1 = Grid3d({ nx, ny, nz }, { dx, dy, dz }, {});

    double val = 0.0;
    for (auto& v : grid1.values())
    {
        v = val;
        val += 1.0;
    }

    const auto grid2 = grid1;

    EXPECT_FLOAT_EQ(
            grid2.at(2, 1, 1),
            grid2.at_pos_pbc(RVec{ 1.25f + (3.0f * box_x), 1.5f - (2.0f * box_y), 3.0f - box_z }, box));
}

} // namespace
} // namespace test
} // namespace gmx
