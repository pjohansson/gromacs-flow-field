/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
/*! \internal \file
 * \brief Tests for H5MD compression routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_compression.h"

#include <hdf5.h>

#include <H5Z_SZ3.hpp>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_datasetbase.h"
#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

constexpr char* enumToString(const H5mdCompression compression)
{
    switch (compression)
    {
        case H5mdCompression::Uncompressed: return "Uncompressed";
        case H5mdCompression::LosslessNoShuffle: return "LosslessNoShuffle";
        case H5mdCompression::LosslessShuffle: return "LosslessShuffle";
        case H5mdCompression::LossySz3BioMD: return "LossySz3BioMD";
        default: GMX_THROW(InvalidInputError("Unknown compression scheme"));
    }
}

//! \brief Helper function for GTest to print test parameters.
void PrintTo(const H5mdCompression& compression, std::ostream* os)
{
    *os << enumToString(compression);
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfTest(const ::testing::TestParamInfo<H5mdCompression>& info)
{
    std::string testName(enumToString(info.param));

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");

    return testName;
}

namespace test
{
namespace
{

template<typename ValueType>
class H5mdSz3CompressionTest : public H5mdTestBase
{
};

// We only compile lossy SZ3 compression for coordinate types.
using TypesToTestForSz3 = ::testing::Types<BasicVector<float>, BasicVector<double>>;
TYPED_TEST_SUITE(H5mdSz3CompressionTest, TypesToTestForSz3);

TYPED_TEST(H5mdSz3CompressionTest, FilterIsApplied)
{
    const auto [creationPropertyList, creationPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_CREATE));

    setSz3Compression<TypeParam>(creationPropertyList, 5e-4);
    EXPECT_EQ(H5Pget_nfilters(creationPropertyList), 1);
    // H5Pget_filter_by_id2 returns: >=0 if the filter (second argument) is set, else <0
    EXPECT_GE(H5Pget_filter_by_id2(
                      creationPropertyList, H5Z_FILTER_SZ3, nullptr, nullptr, nullptr, 0, nullptr, nullptr),
              0);
}

TYPED_TEST(H5mdSz3CompressionTest, ThrowsForNonPositiveAbsoluteError)
{
    const auto [creationPropertyList, creationPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_CREATE));

    EXPECT_THROW(setSz3Compression<TypeParam>(creationPropertyList, -5e-4), H5mdError);
    EXPECT_THROW(setSz3Compression<TypeParam>(creationPropertyList, 0.0), H5mdError);
    EXPECT_NO_THROW(setSz3Compression<TypeParam>(creationPropertyList, +5e-4));
}

TYPED_TEST(H5mdSz3CompressionTest, ThrowsForInvalidHandle)
{
    EXPECT_THROW(setSz3Compression<TypeParam>(H5I_INVALID_HID, 5e-4), H5mdError);
}

//! \brief Test fixture for isCompressedWithSz3() tests
class H5mdIsCompressedWithSz3Test : public H5mdTestBase, public ::testing::WithParamInterface<H5mdCompression>
{
};

TEST_P(H5mdIsCompressedWithSz3Test, Float)
{
    H5mdDataSetBase<BasicVector<float>> dataSet =
            H5mdDataSetBuilder<BasicVector<float>>(this->fileid(), "testDataSet")
                    .withDimension({ 0 })
                    .withCompression(GetParam(), 5e-4)
                    .build();

    if (GetParam() == H5mdCompression::LossySz3BioMD)
    {
        EXPECT_TRUE(isCompressedWithSz3(dataSet.id()));
    }
    else
    {
        EXPECT_FALSE(isCompressedWithSz3(dataSet.id()));
    }
}

TEST_P(H5mdIsCompressedWithSz3Test, Double)
{
    H5mdDataSetBase<BasicVector<double>> dataSet =
            H5mdDataSetBuilder<BasicVector<double>>(this->fileid(), "testDataSet")
                    .withDimension({ 3 })
                    .withCompression(GetParam(), 5e-4)
                    .build();

    if (GetParam() == H5mdCompression::LossySz3BioMD)
    {
        EXPECT_TRUE(isCompressedWithSz3(dataSet.id()));
    }
    else
    {
        EXPECT_FALSE(isCompressedWithSz3(dataSet.id()));
    }
}

INSTANTIATE_TEST_SUITE_P(Works,
                         H5mdIsCompressedWithSz3Test,
                         ::testing::Values(H5mdCompression::Uncompressed,
                                           H5mdCompression::LosslessNoShuffle,
                                           H5mdCompression::LosslessShuffle,
                                           H5mdCompression::LossySz3BioMD),
                         nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
