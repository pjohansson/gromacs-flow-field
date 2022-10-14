/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief Implements position restraint tests.
 *
 * \author Kevin Boyd <kevin44boyd@gmail.com>
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "gromacs/listed_forces/position_restraints.h"

#include <cmath>

#include <memory>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringstream.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! Tolerance for float evaluation
constexpr float c_precisionTolerance = 1e-6;


class RefCoordScalingTest : public ::testing::TestWithParam<std::tuple<RefCoordScaling, PbcType>>
{
public:
    std::vector<RVec> x_;
    std::vector<RVec> f_;

    matrix  box_;
    t_pbc   pbc_;
    PbcType pbcType_;

    t_nrnb                           nrnb_;
    t_forcerec                       fr_;
    InteractionDefinitions           idef_;
    gmx_enerdata_t                   enerd_;
    std::unique_ptr<ForceWithVirial> forceWithVirial_;
    RefCoordScaling                  refCoordScaling_;

    RefCoordScalingTest() : idef_({}), enerd_(1, 0)
    {
        refCoordScaling_ = std::get<0>(GetParam());
        pbcType_         = std::get<1>(GetParam());

        clear_mat(box_);
        box_[0][0] = 0.9;
        box_[1][1] = 1.0;
        box_[2][2] = 1.1;
        set_pbc(&pbc_, pbcType_, box_);

        // FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(1.0, c_precisionTolerance);
        // checker_.setDefaultTolerance(tolerance);

        fr_.rc_scaling    = refCoordScaling_;
        fr_.pbcType       = pbcType_;
    }

    //! Prepares the test with the coordinate and force constant input.
    void setValues(gmx::ArrayRef<const RVec> positions,
                   gmx::ArrayRef<const RVec> referencePositions,
                   gmx::ArrayRef<const RVec> forceConstants,
                   const std::vector<RVec>&  posres_com)
    {
    SCOPED_TRACE(formatString("Setting values!"));
        x_.resize(positions.size());
        std::copy(positions.begin(), positions.end(), x_.begin());

        for (size_t i = 0; i < positions.size(); i++)
        {
            // First item is "type" - each atom will have a different forceparam type
            // Second item is index - we'll just go from 0.
            idef_.il[F_POSRES].iatoms.push_back(i);
            idef_.il[F_POSRES].iatoms.push_back(i);

            auto& entry = idef_.iparams_posres.emplace_back();
            copy_rvec(referencePositions[i], entry.posres.pos0A);
            copy_rvec(forceConstants[i], entry.posres.fcA);
            clear_rvec(entry.posres.pos0B);
            clear_rvec(entry.posres.fcB);
        }
        f_.resize(x_.size(), { 0, 0, 0 });
        forceWithVirial_ = std::make_unique<ForceWithVirial>(f_, /*computeVirial=*/true);

        fr_.posres_com = posres_com;
        fr_.posres_comB = posres_com;
    }
};

// Test setup which mimics the test in position_restraints.cpp for refcoord-scaling = None, All
// Difference is that 
class RefCoordScalingNoneAllTest : public RefCoordScalingTest 
{
public:
    InteractionDefinitions           idef_;
    gmx_enerdata_t                   enerd_;
    TestReferenceData    refData_;
    TestReferenceChecker checker_;

    RefCoordScalingNoneAllTest() : idef_({}), enerd_(1, 0), checker_(refData_.rootChecker())
    {
        refCoordScaling_ = std::get<0>(GetParam());
        pbcType_         = std::get<1>(GetParam());

        clear_mat(box_);
        box_[0][0] = 0.9;
        box_[1][1] = 1.0;
        box_[2][2] = 1.1;
        set_pbc(&pbc_, pbcType_, box_);

        FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(1.0, c_precisionTolerance);
        checker_.setDefaultTolerance(tolerance);

        fr_.rc_scaling    = refCoordScaling_;
        fr_.pbcType       = pbcType_;
    }

    //! Prepares the test with the coordinate and force constant input.
    void setValues(gmx::ArrayRef<const RVec> positions,
                   gmx::ArrayRef<const RVec> referencePositions,
                   gmx::ArrayRef<const RVec> forceConstants,
                   const std::vector<RVec>&  posres_com)
    {
    SCOPED_TRACE(formatString("Setting values!"));
        x_.resize(positions.size());
        std::copy(positions.begin(), positions.end(), x_.begin());

        for (size_t i = 0; i < positions.size(); i++)
        {
            // First item is "type" - each atom will have a different forceparam type
            // Second item is index - we'll just go from 0.
            idef_.il[F_POSRES].iatoms.push_back(i);
            idef_.il[F_POSRES].iatoms.push_back(i);

            auto& entry = idef_.iparams_posres.emplace_back();
            copy_rvec(referencePositions[i], entry.posres.pos0A);
            copy_rvec(forceConstants[i], entry.posres.fcA);
            clear_rvec(entry.posres.pos0B);
            clear_rvec(entry.posres.fcB);
        }
        f_.resize(x_.size(), { 0, 0, 0 });
        forceWithVirial_ = std::make_unique<ForceWithVirial>(f_, /*computeVirial=*/true);

        fr_.posres_com = posres_com;
        fr_.posres_comB = posres_com;
    }
};

std::array<real, static_cast<size_t>(FreeEnergyPerturbationCouplingType::Count)> c_emptyLambdas = { { 0 } };

TEST_P(RefCoordScalingNoneAllTest, ModeNoneOrAllDoesNotUseComGroupInds)
{
    SCOPED_TRACE(formatString("Testing PBC type: %s, refcoord type: %s",
                              c_pbcTypeNames[pbcType_].c_str(),
                              enumValueToString(refCoordScaling_)));
    const std::vector<RVec> positions          = { { 0.0, 0.0, 0.0 }, { 0.4, 0.5, 0.6 } };
    const std::vector<RVec> referencePositions = { { 0.0, 0.0, 0.0 }, { 0.5, 0.6, 0.0 } };
    const std::vector<RVec> forceConstants     = { { 1000, 500, 250 }, { 0, 200, 400 } };

    const std::vector<RVec> posres_com {{0.0, 0.5, 0.0}};
    // these group inds do not exist in posres_com but since they are not used
    // for refcoord-scaling = No,All the test will pass
    const std::vector<unsigned short> refScaleComInds { 100, 100 };

    setValues(positions, referencePositions, forceConstants, posres_com);
    posres_wrapper(&nrnb_,
                   idef_,
                   &pbc_,
                   as_rvec_array(x_.data()),
                   &enerd_,
                   c_emptyLambdas,
                   &fr_,
                   refScaleComInds,
                   numPbcDimensions(pbcType_),
                   forceWithVirial_.get());
    checker_.checkSequence(
            std::begin(forceWithVirial_->force_), std::end(forceWithVirial_->force_), "Forces");
    checker_.checkSequenceArray(3, forceWithVirial_->getVirial(), "Virial contribution");
    checker_.checkReal(enerd_.term[F_POSRES], "Potential energy");
}

TEST_P(RefCoordScalingTest, ComModeUsesComGroupForScaling)
{
    SCOPED_TRACE(formatString("Testing PBC type: %s, refcoord type: %s",
                              c_pbcTypeNames[pbcType_].c_str(),
                              enumValueToString(refCoordScaling_)));
    const std::vector<RVec> positions          = { { 0.1, 0.2, 0.3 }, { 0.4, 0.5, 0.6 } };
    const std::vector<RVec> referencePositions = { { 0.1, 0.0, 0.4 }, { 0.5, 0.6, 0.7 } };
    const std::vector<RVec> forceConstants     = { { 1000, 500, 250 }, { 0, 200, 400 } };

    const std::vector<RVec> posres_com {{0.0, 0.5, 0.0}, {0.0, 1.0, 0.0}};
    // reverse order for atom COM groups to check that the correct one is selected
    const std::vector<unsigned short> refScaleComInds { 1, 0 };

    setValues(positions, referencePositions, forceConstants, posres_com);

    const auto nBoundedDim = numPbcDimensions(pbcType_);
    posres_wrapper(&nrnb_,
                   idef_,
                   &pbc_,
                   as_rvec_array(x_.data()),
                   &enerd_,
                   c_emptyLambdas,
                   &fr_,
                   refScaleComInds,
                   nBoundedDim,
                   forceWithVirial_.get());

    for (size_t i = 0; i < positions.size(); ++i)
    {
        const auto& xs = positions.at(i);
        const auto& ref = referencePositions.at(i);
        const auto& ks = forceConstants.at(i);
        const auto& coms = posres_com.at(refScaleComInds.at(i));

        const auto& fs = forceWithVirial_->force_.at(i);

        for (int d = 0; d < DIM; ++d)
        {
            const auto com_sc = coms[d] * box_[d][d];
            float ref_com;

            if (d < nBoundedDim)
            {
                ref_com = fmod(ref[d] + com_sc, box_[d][d]);
            }
            else
            {
                ref_com = ref[d];
            }

            const auto dx = xs[d] - ref_com;

            const auto force = -ks[d] * dx;
            EXPECT_FLOAT_EQ(fs[d], force);
        }
    }
}

TEST_P(RefCoordScalingTest, RefScaleComGroupIndsEmptyDefaultsToGroup0)
{
    SCOPED_TRACE(formatString("Testing PBC type: %s, refcoord type: %s",
                              c_pbcTypeNames[pbcType_].c_str(),
                              enumValueToString(refCoordScaling_)));
    const std::vector<RVec> positions          = { { 0.1, 0.2, 0.3 }, { 0.4, 0.5, 0.6 } };
    const std::vector<RVec> referencePositions = { { 0.1, 0.0, 0.4 }, { 0.5, 0.6, 0.7 } };
    const std::vector<RVec> forceConstants     = { { 1000, 500, 250 }, { 0, 200, 400 } };

    const std::vector<RVec> posres_com {{0.0, 0.5, 0.0}};
    const std::vector<unsigned short> refScaleComInds;

    setValues(positions, referencePositions, forceConstants, posres_com);

    const auto nBoundedDim = numPbcDimensions(pbcType_);
    posres_wrapper(&nrnb_,
                   idef_,
                   &pbc_,
                   as_rvec_array(x_.data()),
                   &enerd_,
                   c_emptyLambdas,
                   &fr_,
                   refScaleComInds,
                   nBoundedDim,
                   forceWithVirial_.get());

    for (size_t i = 0; i < positions.size(); ++i)
    {
        const auto& xs = positions.at(i);
        const auto& ref = referencePositions.at(i);
        const auto& ks = forceConstants.at(i);
        const auto& coms = posres_com.at(0);

        const auto& fs = forceWithVirial_->force_.at(i);

        for (int d = 0; d < DIM; ++d)
        {
            const auto com_sc = coms[d] * box_[d][d];
            float ref_com;

            if (d < nBoundedDim)
            {
                ref_com = fmod(ref[d] + com_sc, box_[d][d]);
            }
            else
            {
                ref_com = ref[d];
            }

            const auto dx = xs[d] - ref_com;

            const auto force = -ks[d] * dx;
            EXPECT_FLOAT_EQ(fs[d], force);
        }
    }
}

//! PBC values for testing
std::vector<PbcType> c_pbcForTests = { PbcType::No, PbcType::XY, PbcType::Xyz };
//! Reference Coordinate Scaling values for testing
std::vector<RefCoordScaling> c_refCoordScalingModesAllNone = { RefCoordScaling::No,
                                                               RefCoordScaling::All };

INSTANTIATE_TEST_SUITE_P(PosResRefScaleTest,
                         RefCoordScalingNoneAllTest,
                         ::testing::Combine(::testing::ValuesIn(c_refCoordScalingModesAllNone),
                                            ::testing::ValuesIn(c_pbcForTests)));
INSTANTIATE_TEST_SUITE_P(PosResRefScaleTest,
                         RefCoordScalingTest,
                         ::testing::Combine(::testing::ValuesIn(std::vector<RefCoordScaling>{RefCoordScaling::Com}),
                                            ::testing::ValuesIn(c_pbcForTests)));

} // namespace

} // namespace test

} // namespace gmx
