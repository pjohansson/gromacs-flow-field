#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

#include "flow_field.h"
#include "utils.h"

namespace flow
{

/*******************
 * DATA COLLECTION *
 *******************/

static void
add_flow_to_bin(flow::Bin     &bin,
                const rvec     v,
                const real     mass)
{
    bin[FlowVar::NumAtoms] += 1.0;
    bin[FlowVar::Temp    ] += mass * norm2(v);
    bin[FlowVar::Mass    ] += mass;
    bin[FlowVar::U       ] += mass * v[XX];
    bin[FlowVar::V       ] += mass * v[ZZ];
}


static void
collect_flow_data(flow::FlowData         &flowcr,
                  const t_commrec        *cr,
                  const t_inputrec       *ir,
                  const t_mdatoms        *mdatoms,
                  const t_state          *state,
                  const SimulationGroups *groups)
{
    // Atom position buffer
    rvec r;

    // Velocity buffer and length of a half-time-step, to be used
    // if we need to project positions backwards when using the
    // leap-frog integrator
    rvec v_buf;
    const auto dt_half = static_cast<real>(0.5 * ir->delta_t);
    const bool integratorIsLeapFrog = (ir->eI == IntegrationAlgorithm::MD);

    const int num_groups = flowcr.group_data.empty() ? 1 : flowcr.group_data.size();

    for (size_t i = 0; i < static_cast<size_t>(mdatoms->homenr); ++i)
    {
        // Check for match to the input group using the global atom index,
        // since groups contain these indices instead of MPI rank local indices
        const auto index_global = haveDDAtomOrdering(*cr)
            ? cr->dd->globalAtomIndices[i]
            : static_cast<int>(i);

        const auto index_group = getGroupType(*groups, SimulationAtomGroupType::User1, index_global);

        if (index_group < num_groups)
        {
            copy_rvec(state->x[i], r);
            const auto v = state->v[i];
            const auto mass = mdatoms->massT[i];

            /* Fix by Michele Pellegrino */
            /* If we are using the leap-frog integrator, project the positions
               back in time one-half step so that both positions and velocities
               are at the same time. */
            if (integratorIsLeapFrog)
            {
                svmul(dt_half, v, v_buf);
                rvec_dec(r, v_buf);
            }

            auto& bin = flowcr.flow_field.at_pos_pbc(r, state->box);
            add_flow_to_bin(bin, v, mass);

            /* This checks for whether the current atom belongs to a specific
               group, if multiple groups are selected. But, I no longer understand
               exactly what the check does.

               TODO: Figure this out. // Petter */
            if (
                !flowcr.group_data.empty()
                && (index_group < static_cast<int>(flowcr.group_data.size()))
            )
            {
                auto& bin = flowcr.group_data.at(index_group).at_pos_pbc(r, state->box);
                add_flow_to_bin(bin, v, mass);
            }
        }
    }

    ++flowcr.num_samples;
}


/**********
 * OUTPUT *
 **********/

static void
average_flow_field(FlowField &flow_field, const size_t num_samples_int)
{
    const auto num_samples = static_cast<double>(num_samples_int);
    const auto bin_volume = static_cast<double>(
        flow_field.spacing[XX] * flow_field.spacing[YY] * flow_field.spacing[ZZ]
    );

    for (auto& bin : flow_field.values)
    {
        const auto num_atoms = bin[FlowVar::NumAtoms];
        const auto mass      = bin[FlowVar::Mass    ];

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
        bin[FlowVar::Mass]     /= (num_samples * bin_volume);
    }
}

//! A `Bin` dressed with bin indices along x and z in the flow field grid
struct IndexedBin
{
    //! Index along the x axis
    size_t ix;
    //! Index along the z axis
    size_t iz;
    //! Flow field data in bin
    Bin values;
};

//! Data from a flow field which has been prepared for output
struct Output {
    //! Constructor which copies metadata from given `flow_field`
    Output(const FlowField &flow_field)
    :fnbase { flow_field.fnbase },
     shape { flow_field.shape },
     spacing { flow_field.spacing } {}

    //! Base filename for output (`[fnbase]_00001.dat`, ...)
    std::string fnbase;

    //! Grid shape
    gmx::IVec shape;

    //! Grid bin spacing
    gmx::RVec spacing;

    //! Grid origin in system coordinates
    gmx::RVec origin = { 0.0, 0.0, 0.0 };

    //! Non-empty bins with positional indices
    std::vector<IndexedBin> bins;
};

//! Data for all collected flow fields, prepared for output
struct OutputFields {
    //! Main flow field
    Output full;
    //! Sub group flow fields
    std::vector<Output> groups;
};

static Output
get_single_output_flow_field(const FlowField &flow_field)
{
    auto output = Output{flow_field};

    for (size_t ix = 0; ix < flow_field.nx(); ++ix)
    {
        for (size_t iz = 0; iz < flow_field.nz(); ++iz)
        {
            const auto& bin = flow_field.at(ix, 0, iz);

            if (bin[FlowVar::Mass] > 0.0)
            {
                output.bins.push_back(IndexedBin{
                    ix, iz, bin
                });
            }
        }
    }

    return output;
}

static OutputFields
get_averaged_flow_bins(FlowData &flowcr)
{
    average_flow_field(flowcr.flow_field, flowcr.num_samples);
    const auto full_field = get_single_output_flow_field(flowcr.flow_field);

    std::vector<Output> group_fields;
    for (auto& group_field : flowcr.group_data)
    {
        average_flow_field(group_field, flowcr.num_samples);
        group_fields.push_back(get_single_output_flow_field(group_field));
    }

    return OutputFields{full_field, group_fields};
}


static void
write_header(FILE         *fp,
             const size_t  nx,
             const size_t  ny,
             const double  dx,
             const double  dy,
             const size_t  num_values)
{
    std::ostringstream buf;

    buf << "FORMAT " << FLOW_FILE_HEADER_NAME << '\n';
    buf << "ORIGIN 0.0 0.0\n";
    buf << "SHAPE " << nx << ' ' << ny << '\n';
    buf << "SPACING " << dx << ' ' << dy << '\n';
    buf << "NUMDATA " << num_values << '\n';
    buf << "FIELDS IX IY N T M U V\n";
    buf << "COMMENT Grid is regular but only non-empty bins are output\n";
    buf << "COMMENT There are 'NUMDATA' non-empty bins and that many values are stored for each field\n";
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

    const std::string header_str { buf.str() };

    fwrite(header_str.c_str(), sizeof(char), header_str.size(), fp);
}


static void 
write_flow_field_to_disk(const Output &flow_field, const size_t file_index)
{
    char fn[STRLEN];

    snprintf(
        fn, STRLEN,
        "%s_%05lu.%s",
        flow_field.fnbase.c_str(), 
        file_index, 
        ftp2ext(efDAT)
    );

    FILE *fp = gmx_ffopen(fn, "wb");

    write_header(
        fp, 
        flow_field.shape[XX], 
        flow_field.shape[ZZ], 
        flow_field.spacing[XX],
        flow_field.spacing[ZZ],
        flow_field.bins.size()
    );

    std::vector<uint64_t> buf_ix,
                          buf_iz;
    std::vector<float> buf_num_density,
                       buf_mass,
                       buf_temp,
                       buf_vx,
                       buf_vz;

    const auto num_bins = flow_field.bins.size();
    buf_ix.reserve(num_bins);
    buf_iz.reserve(num_bins);
    buf_num_density.reserve(num_bins);
    buf_mass.reserve(num_bins);
    buf_temp.reserve(num_bins);
    buf_vx.reserve(num_bins);
    buf_vz.reserve(num_bins);

    for (const auto& bin : flow_field.bins)
    {
        buf_ix.push_back(bin.ix);
        buf_iz.push_back(bin.iz);
        buf_num_density.push_back(bin.values[FlowVar::NumAtoms]);
        buf_mass.push_back(bin.values[FlowVar::Mass]);
        buf_temp.push_back(bin.values[FlowVar::Temp]);
        buf_vx.push_back(bin.values[FlowVar::U]);
        buf_vz.push_back(bin.values[FlowVar::V]);
    }

    fwrite(buf_ix.data(),           sizeof(uint64_t), num_bins, fp);
    fwrite(buf_iz.data(),           sizeof(uint64_t), num_bins, fp);
    fwrite(buf_num_density.data(),  sizeof(float),    num_bins, fp);
    fwrite(buf_temp.data(),         sizeof(float),    num_bins, fp);
    fwrite(buf_mass.data(),         sizeof(float),    num_bins, fp);
    fwrite(buf_vx.data(),           sizeof(float),    num_bins, fp);
    fwrite(buf_vz.data(),           sizeof(float),    num_bins, fp);

    gmx_ffclose(fp);
}


static void
write_all_flow_fields_to_disk(const OutputFields &output_fields,
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

static void
mpi_collect_flow_data_on_master(FlowData        &flowcr,
                                const t_commrec *cr)
{
    if (PAR(cr))
    {
        MPI_Reduce(MASTER(cr) ? MPI_IN_PLACE : flowcr.flow_field.values.data(),
                MASTER(cr) ? flowcr.flow_field.values.data() : NULL,
                flowcr.flow_field.values.size(),
                MPI_DOUBLE, MPI_SUM, MASTERRANK(cr),
                cr->mpi_comm_mygroup);

        for (auto& group_data : flowcr.group_data)
        {
            MPI_Reduce(MASTER(cr) ? MPI_IN_PLACE : group_data.values.data(),
                    MASTER(cr) ? group_data.values.data() : NULL,
                    group_data.values.size(),
                    MPI_DOUBLE, MPI_SUM, MASTERRANK(cr),
                    cr->mpi_comm_mygroup);
        }
    }
}


/********************
 * PUBLIC FUNCTIONS *
 ********************/

FlowData
init_flow_container(const int               nfile,
                    const t_filenm          fnm[],
                    const t_inputrec       *ir,
                    const SimulationGroups *groups,
                    const t_state          *state)
{
    const auto step_collect = static_cast<uint64_t>(ir->userint1);
    auto step_output = static_cast<uint64_t>(ir->userint2);

    const auto nx = ir->userint3;
    const auto nz = ir->userint4;

    // Control userargs, although this should be done during pre-processing
    if (nx <= 0 || nz <= 0)
    {
        gmx_fatal(FARGS,
                  "Number of bins along x (userint3 = %d) and z (userint4 = %d) "
                  "for flow data calculation and output must be larger than 0.",
                  nx, nz);
    }

    if (step_collect <= 0 || step_output <= 0)
    {
        gmx_fatal(FARGS,
                  "Number of steps that elapse between collection (userint1 = %lu) "
                  "and output (userint2 = %lu) of flow data must be larger than 0.",
                  step_collect, step_output);
    }
    else if (step_collect > step_output)
    {
        gmx_fatal(FARGS,
                  "Number of steps elapsing between output (userint2 = %lu) "
                  "must be larger than steps between collection (userint1 = %lu).",
                  step_output, step_collect);
    }
    else if (step_output % step_collect != 0)
    {
        const auto new_step_output = static_cast<uint64_t>(
            round(step_output / step_collect) * step_collect
        );

        gmx_warning("Steps for outputting flow data (userint2 = %lu) not "
                    "multiple of steps for collecting (userint1 = %lu). "
                    "Setting number of steps that elapse between output to %lu.",
                    step_output, step_collect, new_step_output);

        step_output = new_step_output;
    }

    // Get name base of output datamaps by stripping the extension and dot (.)
    std::string fnbase = opt2fn("-flow", nfile, fnm);

    const int ext_length = static_cast<int>(strlen(ftp2ext(efDAT)));
    const int base_length = static_cast<int>(fnbase.size()) - ext_length - 1;

    if (base_length > 0)
    {
        fnbase.resize(static_cast<size_t>(base_length));
    }

    // If more than one group is selected for output,
    // collect them to do separate collection for each
    // individual group (as well as them all combined)
    const size_t num_groups = get_num_groups(groups);

    std::vector<std::string> group_names;

    if (num_groups > 1)
    {
        for (size_t i = 0; i < num_groups; ++i)
        {
            const auto global_group_index
                = groups->groups[SimulationAtomGroupType::User1].at(i);

            const char *name = *groups->groupNames[global_group_index];
            group_names.push_back(std::string(name));
        }
    }

    return FlowData(
        fnbase, group_names, nx, nz, state->box, step_collect, step_output
    );
}


void
print_flow_collection_information(const FlowData       &flowcr,
                                  const double          dt,
                                  const gmx::MDLogger  &mdlog)
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
        .appendTextFormatted(
            "  Collect: %g ps (every %lu steps).\n",
            flowcr.step_collect * dt, flowcr.step_collect
        )
        .appendTextFormatted(
            "  Output:  %g ps (every %lu steps).",
            flowcr.step_output * dt, flowcr.step_output
        );

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendText("Flow field grid information:\n")
        .appendTextFormatted(
            "  Shape:   %lu x %lu (along x and z)\n",
            flowcr.flow_field.nx(), flowcr.flow_field.nz()
        )
        .appendTextFormatted(
            "  Spacing: %g x %g nm^2",
            flowcr.flow_field.dx(), flowcr.flow_field.dz()
        );

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendTextFormatted(
            "Writing full flow data to files "
            "with base '%s_00001.dat' (...).",
            flowcr.flow_field.fnbase.c_str()
        );

    if (!flowcr.group_data.empty())
    {
        GMX_LOG(mdlog.warning)
            .asParagraph()
            .appendText(
                "Multiple groups selected for flow output. Will collect "
                "individual flow data for each group individually in "
                "addition to the combined field:\n"
            );

        for (const auto& group : flowcr.group_data)
        {
            GMX_LOG(mdlog.warning)
                .appendTextFormatted(
                    "  %s -> '%s_00001.dat' (...)\n",
                    group.name.c_str(),
                    group.fnbase.c_str()
                );
        }
    }

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendText("Have a nice day.");

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendText("****************************************\n")
        .appendText("* END FLOW DATA COLLECTION INFORMATION *\n")
        .appendText("****************************************");
}


void
flow_collect_or_output(FlowData               &flowcr,
                       const uint64_t          current_step,
                       const t_commrec        *cr,
                       const t_inputrec       *ir,
                       const t_mdatoms        *mdatoms,
                       const t_state          *state,
                       const SimulationGroups *groups)
{
    collect_flow_data(flowcr, cr, ir, mdatoms, state, groups);

    if (do_per_step(current_step, flowcr.step_output)
        && (static_cast<int64_t>(current_step) != ir->init_step))
    {
        mpi_collect_flow_data_on_master(flowcr, cr);

        if (MASTER(cr))
        {
            const auto output_data = get_averaged_flow_bins(flowcr);
            write_all_flow_fields_to_disk(
                output_data, current_step, flowcr.step_output
            );
        }

        flowcr.reset_data();
    }
}

} // namespace flow
