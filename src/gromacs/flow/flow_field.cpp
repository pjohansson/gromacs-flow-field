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

namespace flow
{

/*! \brief Get the number of groups in User1

    This is slightly complicated by how Gromacs adds a "rest" group
    to the array of names if the other groups do not add up to all
    atoms in the system. Thus, we detect if the final group is called
    exactly "rest" and if so do not count it as one of the groups. */
static size_t
get_num_groups(const SimulationGroups *groups)
{
    size_t num_groups = 0;

    for (const auto global_group_index
         : groups->groups[SimulationAtomGroupType::User1])
    {
        const auto name = groups->groupNames[global_group_index];

        if (strncmp(*name, "rest", 8) == 0)
        {
            break;
        }

        num_groups++;
    }

    return num_groups;
}


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

    // If more than one group is selected for output, collect them to do separate
    // collection for each individual group (as well as them all combined)
    //
    // Get the number of selected groups, subtract 1 because "rest" is always present
    // const size_t num_groups = groups->grps[egcUser1].nr - 1;
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


struct FlowBinData {
    float mass,
          temp,
          num_atoms,
          u,
          v;
};


struct GroupOutput {
    GroupOutput() = default;

    GroupOutput(const size_t num_elements, const std::string fnbase)
    :fnbase { fnbase }
    {
        ix.reserve(num_elements);
        iy.reserve(num_elements);
        mass_density.reserve(num_elements);
        num_density.reserve(num_elements);
        temp.reserve(num_elements);
        us.reserve(num_elements);
        vs.reserve(num_elements);
    }

    std::vector<uint64_t> ix, iy;
    std::vector<float> mass_density, num_density, temp, us, vs;

    std::string fnbase;
};


struct FlowFieldOutput {
    FlowFieldOutput(const FlowData &flowcr)
    :nx { flowcr.flow_field.nx() },
     nz { flowcr.flow_field.nz() },
     dx { flowcr.flow_field.dx() },
     dz { flowcr.flow_field.dz() },
     all_groups { nx * nz, flowcr.flow_field.fnbase }
    {
        const auto num_bins = nx * nz;

        for (const auto& group_data : flowcr.group_data)
        {
            individual_groups.push_back(
                GroupOutput(num_bins, group_data.fnbase)
            );
        }
    }

    size_t nx, nz;
    double dx, dz;

    GroupOutput all_groups;
    std::vector<GroupOutput> individual_groups;
};


static FlowBinData
calc_values_in_bin(const flow::Bin &bin,
                   const uint64_t   num_samples_int)
{
    const auto num_atoms = bin[FlowVar::NumAtoms];
    const auto mass      = bin[FlowVar::Mass    ];

    /* The temperature and flow is averaged by the sampled number
       of atoms and mass in each bin. To not divide by zero in empty
       bins we take care to check. */
    double flow_x = 0.0,
           flow_z = 0.0,
           temperature = 0.0;

    if (num_atoms > 0.0)
    {
        temperature = bin[FlowVar::Temp] / (2.0 * gmx::c_boltz * num_atoms);
    }

    if (mass > 0.0)
    {
        flow_x = bin[FlowVar::U] / mass;
        flow_z = bin[FlowVar::V] / mass;
    }

    /* In contrast to above, the mass and number of atoms has to be divided by
       the number of samples taken to get their average. */
    const auto num_samples = static_cast<double>(num_samples_int);
    const auto avg_num_atoms = num_atoms / num_samples;
    const auto avg_mass = mass / num_samples;

    FlowBinData bin_data;

    bin_data.num_atoms = static_cast<float>(avg_num_atoms);
    bin_data.mass = static_cast<float>(avg_mass);
    bin_data.temp = static_cast<float>(temperature);
    bin_data.u    = static_cast<float>(flow_x);
    bin_data.v    = static_cast<float>(flow_z);

    return bin_data;
}


static void
add_bin_if_non_empty(GroupOutput       &data,
                     const size_t       ix,
                     const size_t       iy,
                     const double       bin_volume,
                     const FlowBinData &bin_data)
{
    if (bin_data.mass > 0.0)
    {
        data.ix.push_back(static_cast<uint64_t>(ix));
        data.iy.push_back(static_cast<uint64_t>(iy));

        data.mass_density.push_back(bin_data.mass / bin_volume);
        data.num_density.push_back(bin_data.num_atoms / bin_volume);
        data.temp.push_back(bin_data.temp);
        data.us.push_back(bin_data.u);
        data.vs.push_back(bin_data.v);
    }
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
write_flow_data(const GroupOutput &output,
                const size_t       num_file,
                const size_t       nx,
                const size_t       ny,
                const double       dx,
                const double       dy)
{
    char fn[STRLEN];

    snprintf(fn,
             STRLEN,
             "%s_%05lu.%s",
             output.fnbase.c_str(), num_file, ftp2ext(efDAT));

    FILE *fp = gmx_ffopen(fn, "wb");

    const size_t num_to_write = output.ix.size();
    write_header(fp, nx, ny, dx, dy, num_to_write);

    fwrite(output.ix.data(),           sizeof(uint64_t), num_to_write, fp);
    fwrite(output.iy.data(),           sizeof(uint64_t), num_to_write, fp);
    fwrite(output.num_density.data(),  sizeof(float),    num_to_write, fp);
    fwrite(output.temp.data(),         sizeof(float),    num_to_write, fp);
    fwrite(output.mass_density.data(), sizeof(float),    num_to_write, fp);
    fwrite(output.us.data(),           sizeof(float),    num_to_write, fp);
    fwrite(output.vs.data(),           sizeof(float),    num_to_write, fp);

    gmx_ffclose(fp);
}


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


static FlowFieldOutput
get_average_flow_data(FlowData &flowcr)
{
    FlowFieldOutput output(flowcr);

    for (size_t ix = 0; ix < flowcr.flow_field.nx(); ++ix)
    {
        for (size_t iz = 0; iz < flowcr.flow_field.nz(); ++iz)
        {
            const auto bin = flowcr.flow_field.at(ix, 0, iz);
            const auto bin_data = calc_values_in_bin(bin, flowcr.num_samples);

            add_bin_if_non_empty(
                output.all_groups, ix, iz, flowcr.bin_volume, bin_data
            );

            auto group_output = output.individual_groups.begin();
            auto group_data = flowcr.group_data.cbegin();

            while (
                (group_output != output.individual_groups.end())
                && (group_data != flowcr.group_data.cend())
            )
            {
                const auto bin = (*group_data).at(ix, 0, iz);
                const auto group_bin_data = calc_values_in_bin(
                    bin, flowcr.num_samples
                );

                add_bin_if_non_empty(
                    *group_output, ix, iz, flowcr.bin_volume, group_bin_data
                );

                ++group_output;
                ++group_data;
            }
        }
    }

    return output;
}


static void
output_flow_data(const FlowFieldOutput &output,
                 const uint64_t         current_step,
                 const uint64_t         step_output)
{
    const auto file_index = static_cast<size_t>(current_step / step_output);

    write_flow_data(
        output.all_groups, file_index,
        output.nx, output.nz, output.dx, output.dz
    );

    for (const auto& group_data : output.individual_groups)
    {
        write_flow_data(
            group_data, file_index,
            output.nx, output.nz, output.dx, output.dz
        );
    }
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
            const auto output_data = get_average_flow_data(flowcr);
            output_flow_data(output_data, current_step, flowcr.step_output);
        }

        flowcr.reset_data();
    }
}

} // namespace flow
