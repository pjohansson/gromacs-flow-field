#include <array>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"

#ifndef MD_FLOW_FIELD
#define MD_FLOW_FIELD

namespace flow
{

constexpr char FLOW_FILE_HEADER_NAME[] = "GMX_FLOW_2";

// We are using a grid along X and Z so we use a separate enum
// to not confuse our indexing with regular XX, YY and ZZ
enum class GridAxes {
    X,
    Z,
    NumAxes
};
constexpr size_t NUM_FLOW_AXES = static_cast<size_t>(GridAxes::NumAxes);

// Indices for different data in array
enum class FlowVariable {
    NumAtoms,
    Temp,
    Mass,   // Mass in bin (amu)
    U,      // Mass flow along X
    V,      //             and Z
    NumVariables
};
constexpr size_t NUM_FLOW_VARIABLES = static_cast<size_t>(FlowVariable::NumVariables);
using Bin = std::array<double, NUM_FLOW_VARIABLES>;

//! Flow field data and associated metadata
struct GroupFlowData {
    //! Empty constructor
    GroupFlowData() {}

    //! Constructor for full flow field data
    GroupFlowData(const std::string& fnbase,
                  const size_t       num_data)
    :fnbase { fnbase },
     name { "_FULL_" },
     data(num_data, 0.0)
    {}

    //! Constructor which adds `group_name` to `fnbase`
    GroupFlowData(const std::string& fnbase_original,
                  const std::string& group_name,
                  const size_t       num_data)
    :name { group_name },
     data(num_data, 0.0)
    {
        fnbase.append(fnbase_original);
        fnbase.append("_");
        fnbase.append(group_name);
    }

    //! Base for output file names (`fnbase_00001.dat`, ...)
    std::string fnbase;

    //! Name or identifier for group
    std::string name;

    //! Flow field data
    std::vector<double> data;
};

class FlowData {
public:
    //! Whether or not to collect flow field data
    bool bDoFlowCollection = false;

    //! 2D grid data
    GroupFlowData data;
    std::vector<GroupFlowData> group_data; // Similar data for all separate atom groups

    //! Collect flow field data at step multiples of this
    uint64_t step_collect;
    //! Average and output flow field data at step multiples of this
    uint64_t step_output;
    //! # of samples per output (i.e. step_output / step_collect)
    uint64_t step_ratio;

    //! Volume of bins in grid
    double bin_volume;

    //! Empty constructor, turns off flow field collection
    FlowData() {}

    //! Constructor for full flow field
    FlowData(const std::string              &fnbase,
             const std::vector<std::string> &group_names,
             const size_t                    nx,
             const size_t                    nz,
             const double                    dx,
             const double                    dy,
             const double                    dz,
             const uint64_t                  step_collect,
             const uint64_t                  step_output)
    :bDoFlowCollection { true },
     step_collect { step_collect },
     step_output { step_output },
     step_ratio { static_cast<uint64_t>(step_output / step_collect) },
     bin_volume { dx * dy * dz },
     num_bins { nx, nz },
     bin_size { dx, dz },
     inv_bin_size { 1.0 / dx, 1.0 / dz }
     {
        const size_t num_data = nx * nz * NUM_FLOW_VARIABLES;

        data = GroupFlowData(fnbase, num_data);
        for (const auto& name : group_names)
        {
            group_data.push_back(GroupFlowData(fnbase, name, num_data));
        }
     }

    double dx() const { return bin_size[static_cast<size_t>(GridAxes::X)]; }
    double dz() const { return bin_size[static_cast<size_t>(GridAxes::Z)]; }

    double inv_dx() const { return inv_bin_size[static_cast<size_t>(GridAxes::X)]; }
    double inv_dz() const { return inv_bin_size[static_cast<size_t>(GridAxes::Z)]; }

    size_t nx() const { return num_bins[static_cast<size_t>(GridAxes::X)]; }
    size_t nz() const { return num_bins[static_cast<size_t>(GridAxes::Z)]; }

    size_t get_1d_index(const size_t ix, const size_t iz) const
    {
        return (iz * nx() + ix) * NUM_FLOW_VARIABLES;
    }

    size_t get_xbin(const real x) const { return get_bin_from_position(x, nx(), inv_dx()); }
    size_t get_zbin(const real z) const { return get_bin_from_position(z, nz(), inv_dz()); }

    float get_x(const size_t ix) const { return get_position(ix, dx()); }
    float get_z(const size_t iz) const { return get_position(iz, dz()); }

    //! Zero all data for all collected flow fields
    void reset_data() {
        data.data.assign(data.data.size(), 0.0);

        for (auto& group : group_data)
        {
            group.data.assign(group.data.size(), 0.0);
        }
    }

private:
    //! Number of bins per axis
    std::array<size_t, NUM_FLOW_AXES> num_bins;

    //! Bin size per axis
    std::array<double, NUM_FLOW_AXES> bin_size;

    //! Inverted bin size per axis
    std::array<double, NUM_FLOW_AXES> inv_bin_size;

    //! Get bin index from a position along an axis
    size_t get_bin_from_position(const real x, const size_t num_bins, const real inv_bin) const
    {
        auto index = static_cast<int>(floor(x * inv_bin)) % static_cast<int>(num_bins);

        while (index < 0)
        {
            index += num_bins;
        }

        return index;
    }

    //! Get system position of bin with given index along an axis
    float get_position(const size_t index, const float bin_size) const
    {
        return (static_cast<float>(index) + 0.5) * bin_size;
    }
};

//! Prepare and return a container for flow field data
FlowData
init_flow_container(const int               nfile,
                    const t_filenm          fnm[],
                    const t_inputrec       *ir,
                    const SimulationGroups *groups,
                    const t_state          *state);

//! Write information about the flow field collection
void
print_flow_collection_information(const FlowData       &flowcr,
                                  const double          dt,
                                  const gmx::MDLogger  &mdlog);

//! If at a collection or output step, perform actions
void
flow_collect_or_output(FlowData               &flowcr,
                       const uint64_t          step,
                       const t_commrec        *cr,
                       const t_inputrec       *ir,
                       const t_mdatoms        *mdatoms,
                       const t_state          *state,
                       const SimulationGroups *groups);

} // namespace flow

#endif // MD_FLOW_FIELD
