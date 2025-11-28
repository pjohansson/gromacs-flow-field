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

/*! \brief Definitions of flow field module routines
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "flow_field.h"

#include <filesystem>
#include <sstream>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"

#include "utils.h"

namespace gmx
{
namespace flow
{

/***************************
 * FLOWFIELD CLASS METHODS *
 ***************************/

FlowField::FlowField(const std::filesystem::path& basePath, const int numBinsX, const int numBinxZ, const matrix box) :
    basePath_{ basePath },
    name_{ "_ALL_GROUPS_" },
    grid_{ Grid3d<Bin>(IVec{ numBinsX, 1, numBinxZ },
                       RVec{ box[XX][XX] / static_cast<real>(numBinsX),
                             box[YY][YY],
                             box[ZZ][ZZ] / static_cast<real>(numBinxZ) }) }
{
}

FlowField::FlowField(const std::filesystem::path& basePath,
                     const std::string&           groupName,
                     const int                    numBinsX,
                     const int                    numBinsZ,
                     const matrix                 box) :
    name_{ groupName },
    grid_{ Grid3d<Bin>(IVec{ numBinsX, 1, numBinsZ },
                       RVec{ box[XX][XX] / static_cast<real>(numBinsX),
                             box[YY][YY],
                             box[ZZ][ZZ] / static_cast<real>(numBinsZ) }) }
{
    basePath_ += basePath;
    basePath_ += "_";
    basePath_ += groupName;
}

ArrayRef<Bin> FlowField::bins()
{
    return grid_.values();
}

ArrayRef<const Bin> FlowField::bins() const
{
    return grid_.values();
}

void FlowField::updateSimulationBox(const matrix newSimulationBox)
{
    grid_.setBox(newSimulationBox);
}

/**************************
 * FLOWDATA CLASS METHODS *
 **************************/

FlowData::FlowData(const std::filesystem::path&      basePath,
                   const ArrayRef<const std::string> groupNames,
                   const size_t                      numBinsX,
                   const size_t                      numBinsZ,
                   const matrix                      box,
                   const uint64_t                    nstCollect,
                   const uint64_t                    nstOutput) :
    isActive_{ true },
    totalFlowField_{ FlowField(basePath, numBinsX, numBinsZ, box) },
    nstCollect_{ nstCollect },
    nstOutput_{ nstOutput },
    numSamples_{ 0 }
{
    totalFlowField_ = FlowField(basePath, numBinsX, numBinsZ, box);

    for (const std::string& name : groupNames)
    {
        perGroupFlowFields_.push_back(FlowField(basePath, name, numBinsX, numBinsZ, box));
    }
}

FlowField& FlowData::totalFlowField()
{
    return totalFlowField_;
}

const FlowField& FlowData::totalFlowField() const
{
    return totalFlowField_;
}

ArrayRef<FlowField> FlowData::perGroupFlowFields()
{
    return perGroupFlowFields_;
};

ArrayRef<const FlowField> FlowData::perGroupFlowFields() const
{
    return perGroupFlowFields_;
};

void FlowData::reset()
{
    for (Bin& bin : totalFlowField_.bins())
    {
        bin.fill(0.0);
    }

    for (FlowField& groupFlowField : perGroupFlowFields_)
    {
        for (Bin& bin : groupFlowField.bins())
        {
            bin.fill(0.0);
        }
    }

    numSamples_ = 0;
}


/*******************
 * DATA COLLECTION *
 *******************/

//! Add flow contribution from a single atom to a bin
//!
//! Note: Flow along x and z is multiplied by the atom mass here. When
//! averaging the flow inside each bin before writing to disk, we divide
//! by the total mass in the bin. Thus, the flow that we are measuring
//! here is *mass-averaged*.
static void addFlowToBin(Bin& bin, const RVec& velocity, const real mass)
{
    bin[FlowVar::NumAtoms] += 1.0;
    bin[FlowVar::Temp] += mass * norm2(velocity);
    bin[FlowVar::Mass] += mass;
    bin[FlowVar::U] += mass * velocity[XX];
    bin[FlowVar::V] += mass * velocity[ZZ];
}

//! Collect flow field data from all selected atoms in the system
//!
//! The collected variables here are:
//!
//! * Number of atoms in each bin (to be divided by bin volume)
//! * Total atom mass in each bin (to be divided by bin volume)
//! * Temperature in each bin (but here we only add upp the kinetic energy)
//! * Mass flow along x in each bin (to be divided by total mass in bins)
//! * Mass flow along z in each bin (to be divided by total mass in bins)
static void collectFlowData(FlowData*               flowContainer,
                            const t_commrec&        commRec,
                            const t_inputrec&       inputRec,
                            const t_mdatoms&        mdAtoms,
                            const t_state&          state,
                            const SimulationGroups& groups)
{
    // Length of a half-time-step, to be used if we need to project
    // positions backwards when using the leap-frog integrator
    const real dtHalf               = static_cast<real>(0.5 * inputRec.delta_t);
    const bool integratorIsLeapFrog = (inputRec.eI == IntegrationAlgorithm::MD);

    const int numGroups = flowContainer->perGroupFlowFields().empty()
                                  ? 1
                                  : flowContainer->perGroupFlowFields().size();

    flowContainer->totalFlowField().updateSimulationBox(state.box);
    for (FlowField& groupFlowField : flowContainer->perGroupFlowFields())
    {
        groupFlowField.updateSimulationBox(state.box);
    }

    for (int i = 0; i < mdAtoms.homenr; ++i)
    {
        // Check for match to the input group using the global atom index,
        // since groups contain these indices instead of MPI rank local indices
        const int globalAtomIndex = haveDDAtomOrdering(commRec) ? commRec.dd->globalAtomIndices[i]
                                                                : static_cast<int>(i);
        if (const int atomGroupIndexInUser1 =
                    getGroupType(groups, SimulationAtomGroupType::User1, globalAtomIndex);
            atomGroupIndexInUser1 < numGroups)
        {
            RVec        position = state.x[i]; // copy (see leap-frog adjustment below)
            const RVec& velocity = state.v[i]; // reference
            const real  mass     = mdAtoms.massT[i];

            /* If we are using the leap-frog integrator, project the positions
               back in time one-half step so that both positions and velocities
               are at the same time. */
            if (integratorIsLeapFrog)
            {
                position -= dtHalf * velocity;
            }

            const size_t binIndex = flowContainer->totalFlowField().binIndexFromPosition(position);
            {
                Bin& bin = flowContainer->totalFlowField().bins()[binIndex];
                addFlowToBin(bin, velocity, mass);
            }

            // If we are collecting flow field data for multiple groups, we add
            // that here. The `indexGroup` corresponds to the indexing in our
            // collection of flow fields.
            if (atomGroupIndexInUser1 < static_cast<int>(flowContainer->perGroupFlowFields().size()))
            {
                Bin& bin = flowContainer->perGroupFlowFields().at(atomGroupIndexInUser1).bins()[binIndex];
                addFlowToBin(bin, velocity, mass);
            }
        }
    }

    ++flowContainer->numSamples();
}


/*************
 * AVERAGING *
 *************/

//! Average a single flow field bin, in-place
//!
//! Note: This also divides the mass and number of atom fields by the
//! bin volume, making them the mass-and-number densities.
static void averageFlowFieldBin(Bin* bin, const double numSamplesAsDouble, const double binVolume)
{
    const double numAtoms = (*bin)[FlowVar::NumAtoms];
    const double mass     = (*bin)[FlowVar::Mass];

    /* The temperature and flow is averaged by the sampled number
    of atoms and mass in each bin. To not divide by zero in empty
    bins we take care to check. */
    if (numAtoms > 0.0)
    {
        (*bin)[FlowVar::Temp] /= (2.0 * gmx::c_boltz * numAtoms);
    }

    if (mass > 0.0)
    {
        (*bin)[FlowVar::U] /= mass;
        (*bin)[FlowVar::V] /= mass;
    }

    // In contrast to above, the mass and number of atoms has to
    // be divided by the number of samples taken to get their average.
    // Additionally, since we want the mass and atom number densities,
    // divide by the bin volume.
    (*bin)[FlowVar::NumAtoms] /= (numSamplesAsDouble * binVolume);
    (*bin)[FlowVar::Mass] /= (numSamplesAsDouble * binVolume);
}

//! Average the flow field data inside all bins, in-place
//!
//! Note: This also divides the mass and number of atom fields by the
//! bin volume, making them the mass-and-number densities.
static void averageFlowField(FlowField* flowField, const size_t numSamples)
{
    const double numSamplesAsDouble = static_cast<double>(numSamples);
    const double binVolume          = flowField->binVolume();

    for (Bin& bin : flowField->bins())
    {
        averageFlowFieldBin(&bin, numSamplesAsDouble, binVolume);
    }
}


/**********************
 * OUTPUT PREPARATION *
 **********************/

//! Data from a flow field which has been prepared for output
class OutputData
{
public:
    //! Constructor which copies metadata from given `FlowField`, trimming empty bins
    OutputData(const FlowField& flowField) :
        basePath_{ flowField.basePath() }, shape_{ flowField.shape() }, spacing_{ flowField.spacing() }
    {
        const int numBins = shape_[XX] * shape_[YY] * shape_[ZZ];

        ix_.reserve(numBins);
        iz_.reserve(numBins);
        numDensity_.reserve(numBins);
        massDensity_.reserve(numBins);
        temperature_.reserve(numBins);
        ux_.reserve(numBins);
        uz_.reserve(numBins);

        for (int ix = 0; ix < shape_[XX]; ++ix)
        {
            for (int iz = 0; iz < shape_[ZZ]; ++iz)
            {
                const Bin& bin = flowField.binAt(ix, 0, iz);
                if (bin[FlowVar::Mass] > 0.0)
                {
                    addBin(ix, iz, bin);
                }
            }
        }
    }

    //! Return a constant reference to \c basePath_.
    const std::filesystem::path& basePath() const { return basePath_; }
    //! Return a constant reference to \c shape_.
    const IVec& shape() const { return shape_; }
    //! Return a constant reference to \c spacing_.
    const RVec& spacing() const { return spacing_; }
    //! Return a constant reference to \c ix_.
    ArrayRef<const uint64_t> ix() const { return ix_; }
    //! Return a constant reference to \c iz_.
    ArrayRef<const uint64_t> iz() const { return iz_; }
    //! Return a constant reference to \c numDensity_.
    ArrayRef<const float> numDensity() const { return numDensity_; }
    //! Return a constant reference to \c massDensity_.
    ArrayRef<const float> massDensity() const { return massDensity_; }
    //! Return a constant reference to \c temperature_.
    ArrayRef<const float> temperature() const { return temperature_; }
    //! Return a constant reference to \c ux_.
    ArrayRef<const float> ux() const { return ux_; }
    //! Return a constant reference to \c uz_.
    ArrayRef<const float> uz() const { return uz_; }

private:
    //! Add indices and values for a grid bin to the output collections
    void addBin(const uint64_t ixForBin, const uint64_t izForBin, const Bin& bin)
    {
        ix_.push_back(ixForBin);
        iz_.push_back(izForBin);
        numDensity_.push_back(bin[FlowVar::NumAtoms]);
        massDensity_.push_back(bin[FlowVar::Mass]);
        temperature_.push_back(bin[FlowVar::Temp]);
        ux_.push_back(bin[FlowVar::U]);
        uz_.push_back(bin[FlowVar::V]);
    }

    //! Base filename for output (`<basePath>_00001.dat`, ...)
    std::filesystem::path basePath_;

    //! Grid shape
    IVec shape_;

    //! Grid bin spacing
    RVec spacing_;

    // Flow field data, separated by type for writing full arrays to
    // output files. All these data vectors must have an identical
    // number of values, and correspond 1-to-1 to individual bins
    // in the grid.

    //! Bin indices along x
    std::vector<uint64_t> ix_;
    //! Bin indices along z
    std::vector<uint64_t> iz_;
    //! Atom number densities
    std::vector<float> numDensity_;
    //! Mass densities
    std::vector<float> massDensity_;
    //! Temperatures
    std::vector<float> temperature_;
    //! Mass-averaged velocities along x
    std::vector<float> ux_;
    //! Mass-averaged velocities along z
    std::vector<float> uz_;
};

//! Data for all collected flow fields, prepared for output
struct OutputFields
{
    //! Total collected flow field
    OutputData totalFlowFieldData;
    //! Per-group flow field data (if more than one group is collected for)
    std::vector<OutputData> perGroupFlowFieldData;
};

//! Average all flow fields, trim empty bins and return formatted for output
static OutputFields getAveragedFlowFieldsForOutput(FlowData* flowContainer)
{
    averageFlowField(&flowContainer->totalFlowField(), flowContainer->numSamples());
    const OutputData totalFlowField(flowContainer->totalFlowField());

    std::vector<OutputData> perGroupFlowFields;
    for (FlowField& groupFlowField : flowContainer->perGroupFlowFields())
    {
        averageFlowField(&groupFlowField, flowContainer->numSamples());
        perGroupFlowFields.push_back(OutputData(groupFlowField));
    }

    return OutputFields{ totalFlowField, perGroupFlowFields };
}


/**********
 * OUTPUT *
 **********/

//! Write a header with metadata and information into an opened binary file
//!
//! The information is writted as plaintext into the binary file using
//! a buffer. Good text editors can inspect this header, which has some
//! (albeit poor) documentation of how to read the full file.
//!
//! The header ends with a written NULL (`\0`) value.
static void writeHeader(FILE*        fp,
                        const size_t numBinsX,
                        const size_t numBinsY,
                        const double spacingX,
                        const double spacingY,
                        const size_t numValues)
{
    std::ostringstream buf;

    buf << "FORMAT " << FLOW_FILE_HEADER_NAME << '\n';
    buf << "ORIGIN 0.0 0.0\n";
    buf << "SHAPE " << numBinsX << ' ' << numBinsY << '\n';
    buf << "SPACING " << spacingX << ' ' << spacingY << '\n';
    buf << "NUMDATA " << numValues << '\n';
    buf << "FIELDS IX IY N T M U V\n";
    buf << "COMMENT Grid is regular but only non-empty bins are output\n";
    buf << "COMMENT There are 'NUMDATA' non-empty bins and that many values are stored for each "
           "field\n";
    buf << "COMMENT Origin and spacing is given in units of nm\n";
    buf << "COMMENT 'FIELDS' is the different fields for each bin:\n";
    buf << "COMMENT 'IX' and 'IY' are bin indices along x and y respectively\n";
    buf << "COMMENT 'N' is the average atom number density (1/nm^3)\n";
    buf << "COMMENT 'M' is the average mass density (amu/nm^3)\n";
    buf << "COMMENT 'T' is the temperature (K)\n";
    buf << "COMMENT 'U' and 'V' is the mass-averaged flow along x and y respectively (nm/ps)\n";
    buf << "COMMENT Data is stored as 'NUMDATA' counts for each field in 'FIELDS', in order\n";
    buf << "COMMENT 'IX' and 'IY' are 64-bit unsigned integers\n";
    buf << "COMMENT Other fields are 32-bit floating point numbers\n";
    buf << "COMMENT Data begins after '\\0' character\n";
    buf << "COMMENT Example: with 'NUMDATA' = 4 and 'FIELDS' = 'IX IY N T', "
        << "the data following the '\\0' marker is 4 + 4 64-bit integers "
        << "and then 4 + 4 32-bit floating point numbers\n";
    buf << '\0';

    const std::string headerString{ buf.str() };

    fwrite(headerString.c_str(), sizeof(char), headerString.size(), fp);
}

//! Opens a file with a given index and writes the trimmed and prepared flow field data into it
static void writeFlowFieldOutputDataToDisk(const OutputData& outputData, const int64_t nextFileIndex)
{
    std::filesystem::path outputPath = outputData.basePath();
    outputPath += formatString("_%05lld", static_cast<long long>(nextFileIndex));
    outputPath.replace_extension(ftp2ext(efDAT));

    const size_t numBins = outputData.ix().size();

    FILE* fp = gmx_ffopen(outputPath, "wb");

    writeHeader(fp,
                outputData.shape()[XX],
                outputData.shape()[ZZ],
                outputData.spacing()[XX],
                outputData.spacing()[ZZ],
                numBins);

    // The order of writing these fields is *fixed*!
    //  -> IX, IY, NUM_DENSITY, TEMP, MASS_DENSITY, UX, UZ
    //
    // This corresponds to what is written in the header, although
    // one has to take care of keeping that information up-to-date
    // if anything changes in this code.
    fwrite(outputData.ix().data(), sizeof(uint64_t), numBins, fp);
    fwrite(outputData.iz().data(), sizeof(uint64_t), numBins, fp);
    fwrite(outputData.numDensity().data(), sizeof(float), numBins, fp);
    fwrite(outputData.temperature().data(), sizeof(float), numBins, fp);
    fwrite(outputData.massDensity().data(), sizeof(float), numBins, fp);
    fwrite(outputData.ux().data(), sizeof(float), numBins, fp);
    fwrite(outputData.uz().data(), sizeof(float), numBins, fp);

    gmx_ffclose(fp);
}

//! Write all flow fields (full system + groups) to disk for the current step
//!
//! The current step changes which file index will be used for the filename.
//! We divide it by the output frequency to get the index. This accounts for
//! restarts from checkpoints, which retains the step counter from the previous
//! simulation.
static void writeAllFlowFieldToDisk(const OutputFields& outputFields,
                                    const int64_t       currentStep,
                                    const uint64_t      nstOutput)
{
    const int64_t nextFileIndex = currentStep / static_cast<int64_t>(nstOutput);

    writeFlowFieldOutputDataToDisk(outputFields.totalFlowFieldData, nextFileIndex);
    for (const auto& groupFlowField : outputFields.perGroupFlowFieldData)
    {
        writeFlowFieldOutputDataToDisk(groupFlowField, nextFileIndex);
    }
}


/*********************
 * MPI COMMUNICATION *
 *********************/

//! Reduce the data from a single flow field from all ranks to main
static void mpiCollectSingleFlowField(FlowField* flowField, const t_commrec& commRec)
{
    // We want to transmit all `flow::Bin`s in a single MPI communication.
    // The bins are stored in a std::vector, which guarantuees that adjacent
    // members are contiguous in memory.
    //
    // Thus, if we know the exact size of each `Bin` (an array of doubles),
    // and how many of them there are in each vector (identical for every rank),
    // we can easily calculate how many doubles to transmit as:
    //
    // # of doubles = bin.size() * vector.size()
    //
    // However, to transmit them all in one go we must assume that each `Bin`
    // contains only the values themselves, no additional data: it must be
    // exactly sizeof(double) * bin.size() large in memory. std::array does
    // not seem to guarantuee this but tends to be a light wrapper around
    // raw memory, in which case the assumption holds.
    //
    // Thus, we should indeed be able to directly transmit the desired number
    // of doubles from each vector's data storage to the main rank.
    //
    // To be completely sure, we here make a quick check that this assumption
    // is valid for the current compiled program.
    GMX_RELEASE_ASSERT(sizeof(double) * FlowVar::NumVars == sizeof(Bin),
                       "std::vector<flow::Bin> is not contiguous, MPI_Reduce will fail");

    MPI_Reduce(MAIN(&commRec) ? MPI_IN_PLACE : flowField->bins().data(),
               MAIN(&commRec) ? flowField->bins().data() : nullptr,
               FlowVar::NumVars * flowField->bins().size(), // total number of doubles stored in vector
               MPI_DOUBLE,
               MPI_SUM,
               MAINRANK(commRec),
               commRec.mpi_comm_mygroup);
}

//! If we are using MPI, collect all flow field data to the main rank
static void mpiCollectFlowDataOnMain(FlowData* flowContainer, const t_commrec& commRec)
{
    if (PAR(&commRec))
    {
        mpiCollectSingleFlowField(&flowContainer->totalFlowField(), commRec);

        for (FlowField& perGroupFlowField : flowContainer->perGroupFlowFields())
        {
            mpiCollectSingleFlowField(&perGroupFlowField, commRec);
        }
    }
}


/********************
 * PUBLIC FUNCTIONS *
 ********************/

FlowData initFlowContainer(const int               nFile,
                           const t_filenm          fnm[],
                           const t_inputrec&       inputRec,
                           const SimulationGroups& groups,
                           const t_state&          state)
{
    const std::filesystem::path basePath =
            std::filesystem::path(opt2fn("-flow", nFile, fnm)).replace_extension("");

    // If more than one group is selected for output,
    // collect them to do separate collection for each
    // individual group (as well as them all combined)
    const std::vector<std::string> groupNames = getGroupsInUser1(groups);

    return FlowData(basePath,
                    groupNames.size() > 1 ? groupNames : std::vector<std::string>{},
                    static_cast<size_t>(inputRec.flowFieldOptions.numBinsX),
                    static_cast<size_t>(inputRec.flowFieldOptions.numBinsZ),
                    state.box,
                    static_cast<uint64_t>(inputRec.flowFieldOptions.nstSample),
                    static_cast<uint64_t>(inputRec.flowFieldOptions.nstOutput));
}

void printFlowCollectionInformation(const FlowData& flowContainer, const double dt, const MDLogger& mdLog)
{
    // Log to warning level, which prints both to md.log and stdout
    // (info level only writes to md.log)

    GMX_LOG(mdLog.warning)
            .asParagraph()
            .appendText("************************************\n")
            .appendText("* FLOW DATA COLLECTION INFORMATION *\n")
            .appendText("************************************");

    GMX_LOG(mdLog.warning)
            .asParagraph()
            .appendText("Flow field collection frequency:\n")
            .appendTextFormatted("  Collect: %g ps (every %llu steps).\n",
                                 flowContainer.nstCollect() * dt,
                                 static_cast<unsigned long long>(flowContainer.nstCollect()))
            .appendTextFormatted("  Output:  %g ps (every %llu steps).",
                                 flowContainer.nstOutput() * dt,
                                 static_cast<unsigned long long>(flowContainer.nstOutput()));

    GMX_LOG(mdLog.warning)
            .asParagraph()
            .appendText("Flow field grid information:\n")
            .appendTextFormatted("  Shape:   %d x %d (along x and z)\n",
                                 flowContainer.totalFlowField().shape()[XX],
                                 flowContainer.totalFlowField().shape()[ZZ])
            .appendTextFormatted("  Spacing: %g x %g nm^2",
                                 flowContainer.totalFlowField().spacing()[XX],
                                 flowContainer.totalFlowField().spacing()[ZZ]);

    GMX_LOG(mdLog.warning)
            .asParagraph()
            .appendTextFormatted(
                    "Writing full flow data to files "
                    "with base '%s_00001.dat' (...).",
                    flowContainer.totalFlowField().basePath().c_str());

    if (!flowContainer.perGroupFlowFields().empty())
    {
        GMX_LOG(mdLog.warning)
                .asParagraph()
                .appendText(
                        "Multiple groups selected for flow output. Will collect "
                        "individual flow data for each group individually in "
                        "addition to the combined field:\n");

        for (const auto& group : flowContainer.perGroupFlowFields())
        {
            GMX_LOG(mdLog.warning)
                    .appendTextFormatted("  %s -> '%s_00001.dat' (...)\n",
                                         group.groupName(),
                                         group.basePath().c_str());
        }
    }

    GMX_LOG(mdLog.warning).asParagraph().appendText("Have a nice day.");

    GMX_LOG(mdLog.warning)
            .asParagraph()
            .appendText("****************************************\n")
            .appendText("* END FLOW DATA COLLECTION INFORMATION *\n")
            .appendText("****************************************");
}

void collectOrOutputFlowFieldData(FlowData*               flowContainer,
                                  const int64_t           currentStep,
                                  const t_commrec&        commRec,
                                  const t_inputrec&       inputRec,
                                  const t_mdatoms&        mdAtoms,
                                  const t_state&          state,
                                  const SimulationGroups& groups,
                                  gmx_wallcycle*          wallCycleCounters)
{
    wallcycle_start(wallCycleCounters, WallCycleCounter::FlowField);

    wallcycle_sub_start(wallCycleCounters, WallCycleSubCounter::FlowFieldCollect);
    collectFlowData(flowContainer, commRec, inputRec, mdAtoms, state, groups);
    wallcycle_sub_stop(wallCycleCounters, WallCycleSubCounter::FlowFieldCollect);

    if (do_per_step(currentStep, flowContainer->nstOutput()) && (currentStep != inputRec.init_step))
    {
        wallcycle_sub_start(wallCycleCounters, WallCycleSubCounter::FlowFieldOutput);
        mpiCollectFlowDataOnMain(flowContainer, commRec);

        if (MAIN(&commRec))
        {
            const OutputFields outputData = getAveragedFlowFieldsForOutput(flowContainer);
            writeAllFlowFieldToDisk(outputData, currentStep, flowContainer->nstOutput());
        }

        flowContainer->reset();
        wallcycle_sub_stop(wallCycleCounters, WallCycleSubCounter::FlowFieldOutput);
    }

    wallcycle_stop(wallCycleCounters, WallCycleCounter::FlowField);
}

} // namespace flow
} // namespace gmx
