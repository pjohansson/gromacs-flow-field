#include "flow_field.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sstream>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
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

FlowField::FlowField(const std::string& fnbase, const int nx, const int nz, const matrix box) :
    fnbase{ fnbase }, name{ "_FULL_" }
{
    _setup_grid_and_finalize(nx, nz, box);
}


FlowField::FlowField(const std::string& fnbase_original,
                     const std::string& group_name,
                     const int          nx,
                     const int          nz,
                     const matrix       box) :
    name{ group_name }
{
    fnbase.append(fnbase_original);
    fnbase.append("_");
    fnbase.append(group_name);

    _setup_grid_and_finalize(nx, nz, box);
}


size_t FlowField::index_from_pos_2d(const real x, const real z) const
{
    auto ix = static_cast<int>(floor(x * invSpacing()[XX])) % shape()[XX];
    while (ix < 0)
    {
        ix += shape()[XX];
    }

    auto iz = static_cast<int>(floor(z * invSpacing()[ZZ])) % shape()[ZZ];
    while (iz < 0)
    {
        iz += shape()[ZZ];
    }

    // grid is zyx ordered and iy = 0, ny = 1, so:
    // iz + iy * nz + ix * (ny * nz) = iz + ix * nz
    return static_cast<size_t>(iz + (ix * shape()[ZZ]));
}


void FlowField::_setup_grid_and_finalize(const int nx, const int nz, const matrix box)
{
    // shape_ = gmx::IVec{ nx, 1, nz };

    // spacing_ = gmx::RVec{ box[XX][XX] / static_cast<real>(nx),
    //                      box[YY][YY],
    //                      box[ZZ][ZZ] / static_cast<real>(nz) };

    _finalize();
}


/**************************
 * FLOWDATA CLASS METHODS *
 **************************/

FlowData::FlowData(const std::string&              fnbase,
                   const std::vector<std::string>& group_names,
                   const size_t                    nx,
                   const size_t                    nz,
                   const matrix                    box,
                   const uint64_t                  step_collect,
                   const uint64_t                  step_output) :
    bDoFlowCollection{ true },
    flow_field{ FlowField(fnbase, nx, nz, box) },
    step_collect{ step_collect },
    step_output{ step_output },
    num_samples{ 0 }
{
    // flow_field = FlowField(fnbase, nx, nz, box);

    for (const auto& name : group_names)
    {
        group_data.push_back(FlowField(fnbase, name, nx, nz, box));
    }
}


void FlowData::reset_data()
{
    for (auto& bin : flow_field.values())
    {
        bin.fill(0.0);
    }

    for (auto& group : group_data)
    {
        for (auto& bin : group.values())
        {
            bin.fill(0.0);
        }
    }

    num_samples = 0;
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
static void add_flow_to_bin(flow::Bin& bin, const rvec v, const real mass)
{
    bin[FlowVar::NumAtoms] += 1.0;
    bin[FlowVar::Temp] += mass * norm2(v);
    bin[FlowVar::Mass] += mass;
    bin[FlowVar::U] += mass * v[XX];
    bin[FlowVar::V] += mass * v[ZZ];
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
static void collect_flow_data(flow::FlowData&         flowcr,
                              const t_commrec*        cr,
                              const t_inputrec*       ir,
                              const t_mdatoms*        mdatoms,
                              const t_state*          state,
                              const SimulationGroups* groups)
{
    // Atom position buffer
    rvec r;

    // Length of a half-time-step, to be used if we need to project
    // positions backwards when using the leap-frog integrator
    const auto dt_half              = static_cast<real>(0.5 * ir->delta_t);
    const bool integratorIsLeapFrog = (ir->eI == IntegrationAlgorithm::MD);

    const int num_groups = flowcr.group_data.empty() ? 1 : flowcr.group_data.size();

    for (size_t i = 0; i < static_cast<size_t>(mdatoms->homenr); ++i)
    {
        // Check for match to the input group using the global atom index,
        // since groups contain these indices instead of MPI rank local indices
        const auto index_global =
                haveDDAtomOrdering(*cr) ? cr->dd->globalAtomIndices[i] : static_cast<int>(i);

        const auto index_group = getGroupType(*groups, SimulationAtomGroupType::User1, index_global);

        if (index_group < num_groups)
        {
            r[XX] = state->x[i][XX];
            r[ZZ] = state->x[i][ZZ];

            const auto v    = state->v[i];
            const auto mass = mdatoms->massT[i];

            /* Fix by Michele Pellegrino */
            /* If we are using the leap-frog integrator, project the positions
               back in time one-half step so that both positions and velocities
               are at the same time. */
            if (integratorIsLeapFrog)
            {
                r[XX] -= dt_half * v[XX];
                r[ZZ] -= dt_half * v[ZZ];
            }

            const size_t bin_index = flowcr.flow_field.index_from_pos_2d(r[XX], r[ZZ]);

            auto& bin = flowcr.flow_field.values().at(bin_index);
            add_flow_to_bin(bin, v, mass);

            // If we are collecting flow field data for multiple groups, we add
            // that here. The `index_group` corresponds to the indexing in our
            // collection of flow fields.
            if (index_group < static_cast<int>(flowcr.group_data.size()))
            {
                auto& bin = flowcr.group_data.at(index_group).values().at(bin_index);
                add_flow_to_bin(bin, v, mass);
            }
        }
    }

    ++flowcr.num_samples;
}


/*************
 * AVERAGING *
 *************/

//! Average a single flow field bin, in-place
//!
//! Note: This also divides the mass and number of atom fields by the
//! bin volume, making them the mass-and-number densities.
static void average_flow_field_bin(Bin& bin, const double num_samples, const double bin_volume)
{
    const auto num_atoms = bin[FlowVar::NumAtoms];
    const auto mass      = bin[FlowVar::Mass];

    /* The temperature and flow is averaged by the sampled number
    of atoms and mass in each bin. To not divide by zero in empty
    bins we take care to check. */
    if (num_atoms > 0.0)
    {
        bin[FlowVar::Temp] /= (2.0 * gmx::c_boltz * num_atoms);
    }

    if (mass > 0.0)
    {
        bin[FlowVar::U] /= mass;
        bin[FlowVar::V] /= mass;
    }

    // In contrast to above, the mass and number of atoms has to
    // be divided by the number of samples taken to get their average.
    // Additionally, since we want the mass and atom number densities,
    // divide by the bin volume.
    bin[FlowVar::NumAtoms] /= (num_samples * bin_volume);
    bin[FlowVar::Mass] /= (num_samples * bin_volume);
}

//! Average the flow field data inside all bins, in-place
//!
//! Note: This also divides the mass and number of atom fields by the
//! bin volume, making them the mass-and-number densities.
static void average_flow_field(FlowField& flow_field, const size_t num_samples_int)
{
    const auto num_samples = static_cast<double>(num_samples_int);
    const auto bin_volume  = static_cast<double>(flow_field.bin_volume());

    for (auto& bin : flow_field.values())
    {
        average_flow_field_bin(bin, num_samples, bin_volume);
    }
}


/**********************
 * OUTPUT PREPARATION *
 **********************/

//! Data from a flow field which has been prepared for output
struct Output
{
    //! Constructor which copies metadata from given `flow_field`
    Output(const FlowField& flow_field) :
        fnbase{ flow_field.fnbase }, shape{ flow_field.shape() }, spacing{ flow_field.spacing() }
    {
        const auto num_bins = shape[XX] * shape[YY] * shape[ZZ];

        ix.reserve(num_bins);
        iz.reserve(num_bins);
        num_density.reserve(num_bins);
        mass_density.reserve(num_bins);
        temperature.reserve(num_bins);
        ux.reserve(num_bins);
        uz.reserve(num_bins);
    }

    //! Add indices and values for a grid bin to the output collections
    void add_bin(const size_t ix_index, const size_t iz_index, const Bin& bin)
    {
        ix.push_back(ix_index);
        iz.push_back(iz_index);
        num_density.push_back(bin[FlowVar::NumAtoms]);
        mass_density.push_back(bin[FlowVar::Mass]);
        temperature.push_back(bin[FlowVar::Temp]);
        ux.push_back(bin[FlowVar::U]);
        uz.push_back(bin[FlowVar::V]);
    }

    //! Base filename for output (`[fnbase]_00001.dat`, ...)
    std::string fnbase;

    //! Grid shape
    gmx::IVec shape;

    //! Grid bin spacing
    gmx::RVec spacing;

    //! Grid origin in system coordinates
    gmx::RVec origin = { 0.0, 0.0, 0.0 };


    // Flow field data, separated by type for writing full arrays to
    // output files. All these data vectors must have an identical
    // number of values, and correspond 1-to-1 to individual bins
    // in the grid.

    //! Bin indices along x
    std::vector<size_t> ix;
    //! Bin indices along z
    std::vector<size_t> iz;
    //! Atom number densities
    std::vector<float> num_density;
    //! Mass densities
    std::vector<float> mass_density;
    //! Temperatures
    std::vector<float> temperature;
    //! Mass-averaged velocities along x
    std::vector<float> ux;
    //! Mass-averaged velocities along z
    std::vector<float> uz;
};


//! Data for all collected flow fields, prepared for output
struct OutputFields
{
    //! Main flow field
    Output full;
    //! Sub group flow fields
    std::vector<Output> groups;
};


//! Collect non-empty bins from a flow field and prepare for output
static Output get_single_output_flow_field(const FlowField& flow_field)
{
    auto output = Output{ flow_field };

    for (size_t ix = 0; ix < flow_field.nx(); ++ix)
    {
        for (size_t iz = 0; iz < flow_field.nz(); ++iz)
        {
            const auto& bin = flow_field.at(ix, 0, iz);

            if (bin[FlowVar::Mass] > 0.0)
            {
                output.add_bin(ix, iz, bin);
            }
        }
    }

    return output;
}


//! Average all flow fields, trim empty bins and return formatted for output
static OutputFields get_averaged_flow_fields_for_output(FlowData& flowcr)
{
    average_flow_field(flowcr.flow_field, flowcr.num_samples);
    const auto full_field = get_single_output_flow_field(flowcr.flow_field);

    std::vector<Output> group_fields;
    for (auto& group_field : flowcr.group_data)
    {
        average_flow_field(group_field, flowcr.num_samples);
        group_fields.push_back(get_single_output_flow_field(group_field));
    }

    return OutputFields{ full_field, group_fields };
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
static void write_header(FILE* fp, const size_t nx, const size_t ny, const double dx, const double dy, const size_t num_values)
{
    std::ostringstream buf;

    buf << "FORMAT " << FLOW_FILE_HEADER_NAME << '\n';
    buf << "ORIGIN 0.0 0.0\n";
    buf << "SHAPE " << nx << ' ' << ny << '\n';
    buf << "SPACING " << dx << ' ' << dy << '\n';
    buf << "NUMDATA " << num_values << '\n';
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

    const std::string header_str{ buf.str() };

    fwrite(header_str.c_str(), sizeof(char), header_str.size(), fp);
}


//! Opens a file with a given index and writes the flow field data into it
static void write_flow_field_to_disk(const Output& flow_field, const size_t file_index)
{
    char fn[STRLEN];

    snprintf(fn, STRLEN, "%s_%05lu.%s", flow_field.fnbase.c_str(), file_index, ftp2ext(efDAT));

    const size_t num_bins = flow_field.ix.size();

    // Why not ensure that we are not writing garbage? I'm pretty sure we
    // are not, but checking costs nothing
    const bool allArraySizesAreEqual =
            ((num_bins == flow_field.iz.size()) && (num_bins == flow_field.num_density.size())
             && (num_bins == flow_field.mass_density.size()) && (num_bins == flow_field.temperature.size())
             && (num_bins == flow_field.ux.size()) && (num_bins == flow_field.uz.size()));

    GMX_RELEASE_ASSERT(allArraySizesAreEqual,
                       "[FLOW_FIELD] Not all flow field containers had the same number "
                       "of non-empty bins\n");


    FILE* fp = gmx_ffopen(fn, "wb");

    write_header(
            fp, flow_field.shape[XX], flow_field.shape[ZZ], flow_field.spacing[XX], flow_field.spacing[ZZ], num_bins);

    // The order of writing these fields is *fixed*!
    //  -> IX, IY, NUM_DENSITY, TEMP, MASS_DENSITY, UX, UZ
    //
    // This corresponds to what is written in the header, although
    // one has to take care of keeping that information up-to-date
    // if anything changes in this code.
    fwrite(flow_field.ix.data(), sizeof(uint64_t), num_bins, fp);
    fwrite(flow_field.iz.data(), sizeof(uint64_t), num_bins, fp);
    fwrite(flow_field.num_density.data(), sizeof(float), num_bins, fp);
    fwrite(flow_field.temperature.data(), sizeof(float), num_bins, fp);
    fwrite(flow_field.mass_density.data(), sizeof(float), num_bins, fp);
    fwrite(flow_field.ux.data(), sizeof(float), num_bins, fp);
    fwrite(flow_field.uz.data(), sizeof(float), num_bins, fp);

    gmx_ffclose(fp);
}


//! Write all flow fields (full system + groups) to disk for the current step
//!
//! The current step changes which file index will be used for the filename.
//! We divide it by the output frequency to get the index. This accounts for
//! restarts from checkpoints, which retains the step counter from the previous
//! simulation.
static void write_all_flow_fields_to_disk(const OutputFields& output_fields,
                                          const uint64_t      step,
                                          const uint64_t      step_output)
{
    const auto file_index = static_cast<size_t>(step / step_output);

    write_flow_field_to_disk(output_fields.full, file_index);
    for (const auto& group_field : output_fields.groups)
    {
        write_flow_field_to_disk(group_field, file_index);
    }
}


/*********************
 * MPI COMMUNICATION *
 *********************/

//! Reduce the data from a single flow field from all ranks to main
static void mpi_collect_single_flow_field(FlowField& flow_field, const t_commrec* cr)
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

    MPI_Reduce(MAIN(cr) ? MPI_IN_PLACE : flow_field.values().data(),
               MAIN(cr) ? flow_field.values().data() : nullptr,
               FlowVar::NumVars * flow_field.values().size(), // total number of doubles stored in vector
               MPI_DOUBLE,
               MPI_SUM,
               MAINRANK(cr),
               cr->mpi_comm_mygroup);
}

//! If we are using MPI, collect all flow field data to the main rank
static void mpi_collect_flow_data_on_master(FlowData& flowcr, const t_commrec* cr)
{
    if (PAR(cr))
    {
        mpi_collect_single_flow_field(flowcr.flow_field, cr);

        for (auto& group_data : flowcr.group_data)
        {
            mpi_collect_single_flow_field(group_data, cr);
        }
    }
}


/********************
 * PUBLIC FUNCTIONS *
 ********************/

FlowData init_flow_container(const int               nfile,
                             const t_filenm          fnm[],
                             const t_inputrec*       ir,
                             const SimulationGroups* groups,
                             const t_state*          state)
{
    // Get name base of output datamaps by stripping the extension and dot (.)
    std::string fnbase = opt2fn("-flow", nfile, fnm);

    const int ext_length  = static_cast<int>(strlen(ftp2ext(efDAT)));
    const int base_length = static_cast<int>(fnbase.size()) - ext_length - 1;

    if (base_length > 0)
    {
        fnbase.resize(static_cast<size_t>(base_length));
    }

    // If more than one group is selected for output,
    // collect them to do separate collection for each
    // individual group (as well as them all combined)
    const size_t num_groups = ::flow::get_num_groups(groups);

    std::vector<std::string> group_names;

    if (num_groups > 1)
    {
        for (size_t i = 0; i < num_groups; ++i)
        {
            const auto global_group_index = groups->groups[SimulationAtomGroupType::User1].at(i);

            const char* name = *groups->groupNames[global_group_index];
            group_names.push_back(std::string(name));
        }
    }

    return FlowData(fnbase,
                    group_names,
                    static_cast<size_t>(ir->flowFieldOptions.nx),
                    static_cast<size_t>(ir->flowFieldOptions.nz),
                    state->box,
                    static_cast<uint64_t>(ir->flowFieldOptions.nstsample),
                    static_cast<uint64_t>(ir->flowFieldOptions.nstoutput));
}


void print_flow_collection_information(const FlowData& flowcr, const double dt, const gmx::MDLogger& mdlog)
{
    // Log to warning level, which prints both to md.log and stdout
    // (info level only writes to md.log)

    GMX_LOG(mdlog.warning)
            .asParagraph()
            .appendText("************************************\n")
            .appendText("* FLOW DATA COLLECTION INFORMATION *\n")
            .appendText("************************************");

    GMX_LOG(mdlog.warning)
            .asParagraph()
            .appendText("Flow field collection frequency:\n")
            .appendTextFormatted("  Collect: %g ps (every %llu steps).\n",
                                 flowcr.step_collect * dt,
                                 static_cast<unsigned long long>(flowcr.step_collect))
            .appendTextFormatted("  Output:  %g ps (every %llu steps).",
                                 flowcr.step_output * dt,
                                 static_cast<unsigned long long>(flowcr.step_output));

    GMX_LOG(mdlog.warning)
            .asParagraph()
            .appendText("Flow field grid information:\n")
            .appendTextFormatted("  Shape:   %lu x %lu (along x and z)\n",
                                 flowcr.flow_field.nx(),
                                 flowcr.flow_field.nz())
            .appendTextFormatted(
                    "  Spacing: %g x %g nm^2", flowcr.flow_field.dx(), flowcr.flow_field.dz());

    GMX_LOG(mdlog.warning)
            .asParagraph()
            .appendTextFormatted(
                    "Writing full flow data to files "
                    "with base '%s_00001.dat' (...).",
                    flowcr.flow_field.fnbase.c_str());

    if (!flowcr.group_data.empty())
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "Multiple groups selected for flow output. Will collect "
                        "individual flow data for each group individually in "
                        "addition to the combined field:\n");

        for (const auto& group : flowcr.group_data)
        {
            GMX_LOG(mdlog.warning)
                    .appendTextFormatted(
                            "  %s -> '%s_00001.dat' (...)\n", group.name.c_str(), group.fnbase.c_str());
        }
    }

    GMX_LOG(mdlog.warning).asParagraph().appendText("Have a nice day.");

    GMX_LOG(mdlog.warning)
            .asParagraph()
            .appendText("****************************************\n")
            .appendText("* END FLOW DATA COLLECTION INFORMATION *\n")
            .appendText("****************************************");
}


void flow_collect_or_output(FlowData&               flowcr,
                            const int64_t           current_step,
                            const t_commrec*        cr,
                            const t_inputrec*       ir,
                            const t_mdatoms*        mdatoms,
                            const t_state*          state,
                            const SimulationGroups* groups)
{
    collect_flow_data(flowcr, cr, ir, mdatoms, state, groups);

    if (do_per_step(current_step, flowcr.step_output) && (current_step != ir->init_step))
    {
        mpi_collect_flow_data_on_master(flowcr, cr);

        if (MAIN(cr))
        {
            const auto output_data = get_averaged_flow_fields_for_output(flowcr);

            write_all_flow_fields_to_disk(
                    output_data, static_cast<uint64_t>(current_step), flowcr.step_output);
        }

        flowcr.reset_data();
    }
}

} // namespace flow
} // namespace gmx
