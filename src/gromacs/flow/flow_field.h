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

/*! \brief Declarations of flow field module routines
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FLOW_FIELD_H
#define GMX_FLOW_FIELD_H

#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include "grid.h"

struct t_filenm;
struct t_inputrec;
struct SimulationGroups;
struct t_state;
struct t_commrec;
struct t_mdatoms;
struct gmx_wallcycle;

namespace gmx
{

template<typename T>
class ArrayRef;
struct MDLogger;

namespace flow
{

//! Unique identifier for flow field data files
constexpr char FLOW_FILE_HEADER_NAME[] = "GMX_FLOW_2";


//! Scoped space for an enum of variables to collect in bins
//!
//! Used to index inside each `Bin` array. The namespace is used
//! to scope each enum variable (good!) but avoid having to always
//! static_cast to integer when using it (bad!).
namespace FlowVar
{
enum FlowVariable : int
{
    //! Number of atoms
    NumAtoms,
    //! Temperature
    Temp,
    //! Mass (in atomic mass units)
    Mass,
    //! Mass flow along x
    U,
    //! Mass flow along z
    V,
    //! Total number of variables
    NumVars
};
} // namespace FlowVar

//! Data stored in a single bin of the flow field grid
using Bin = std::array<double, FlowVar::NumVars>;

//! Flow field data and associated metadata
class FlowField
{
public:
    //! Empty constructor
    FlowField() {}

    //! Constructor for full flow field data
    FlowField(const std::filesystem::path& basePath, int numBinsX, int numBinsZ, const matrix box);

    //! Constructor which joins \p groupName with fnBaseName
    FlowField(const std::filesystem::path& basePath,
              const std::string&           groupName,
              int                          numBinsX,
              int                          numBinsZ,
              const matrix                 box);

    //! Get the bin volume.
    double binVolume() const { return grid_.binVolume(); }

    //! Get the grid shape.
    const IVec& shape() const { return grid_.shape(); }

    //! Get the grid spacing.
    const RVec& spacing() const { return grid_.spacing(); }

    //! Get a reference to the 1d array of bins in the grid.
    ArrayRef<Bin>       bins();
    ArrayRef<const Bin> bins() const;

    //! Get a reference to the bin at a position in the grid (throws if indices are outside).
    const Bin& binAt(const size_t ix, const size_t iy, const size_t iz) const
    {
        return grid_.at(ix, iy, iz);
    }

    //! Get the 1d bin index for the given \p positions (PBC corrected to be put inside the box)
    size_t binIndexFromPosition(rvec position) const { return grid_.indexFromPosition(position); }

    //! Return the base path which flow field files will be written as.
    const std::filesystem::path& basePath() const { return basePath_; }

    //! Return the group name.
    const char* groupName() const { return name_.c_str(); }

    //! Update the simulation box size.
    void updateSimulationBox(const matrix newSimulationBox);

private:
    //! Base for output file names (`fnbase_00001.dat`, ...)
    std::filesystem::path basePath_;

    //! Name or identifier for group
    std::string name_;

    //! Grid of collected values.
    Grid3d<Bin> grid_;
};

struct FlowData
{
    //! Whether or not to collect flow field data
    bool bDoFlowCollection = false;

    //! 2d flow field grid data for all groups (combined field)
    FlowField totalFlowField;

    //! 2d flow field data for individual groups, if multiple are selected
    std::vector<FlowField> perGroupFlowFields;

    //! Collect flow field data at step multiples of this
    uint64_t nstCollect;

    //! Average and output flow field data at step multiples of this
    uint64_t nstOutput;

    //! Number of samples since last output of flow field
    uint64_t numSamples;

    //! Empty constructor, turns off flow field collection
    FlowData() {}

    //! Constructor for full flow field, turns on collection
    FlowData(const std::filesystem::path& basePath,
             ArrayRef<const std::string>  groupNames,
             size_t                       numBinsX,
             size_t                       numBinsZ,
             const matrix                 box,
             uint64_t                     nstCollect,
             uint64_t                     nstOutput);

    //! Zero all data for all collected flow fields and reset sample counter
    void reset();
};

//! Prepare and return a container for flow field data
FlowData initFlowContainer(int                     nFile,
                           const t_filenm          fnm[],
                           const t_inputrec&       inputRecord,
                           const SimulationGroups& groups,
                           const t_state&          state);


//! Write information about the flow field collection
void printFlowCollectionInformation(const FlowData& flowContainer, double dt, const gmx::MDLogger& mdLog);


//! If at a collection or output step, perform actions
void collectOrOutputFlowFieldData(FlowData&               flowContainer,
                                  int64_t                 currentStep,
                                  const t_commrec&        commRec,
                                  const t_inputrec&       inputRec,
                                  const t_mdatoms&        mdAtoms,
                                  const t_state&          state,
                                  const SimulationGroups& groups,
                                  gmx_wallcycle*          wallCycleCounters);

} // namespace flow
} // namespace gmx

#endif // GMX_FLOW_FIELD_H
