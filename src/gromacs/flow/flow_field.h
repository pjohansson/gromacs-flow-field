#include <array>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"

#include "grid.h"

#ifndef MD_FLOW_FIELD
#define MD_FLOW_FIELD

namespace flow
{

constexpr char FLOW_FILE_HEADER_NAME[] = "GMX_FLOW_2";

// Indices for different data in array
enum class FlowVariable : size_t {
    NumAtoms,
    Temp,
    Mass,   // Mass in bin (amu)
    U,      // Mass flow along X
    V,      //             and Z
    NumVariables
};
using Bin = std::array<double, static_cast<size_t>(FlowVariable::NumVariables)>;

//! Flow field data and associated metadata
//!
//! We subclass `Grid3d` in order to deal with the grid bookkeeping.
//! Since the flow field data is more or less an extension of that
//! class we do this rather than creating a `Grid3d` member variable,
//! which makes working with this class more painful and harder to
//! understand.
class GroupFlowData : public Grid3d<Bin> {
public:
    //! Empty constructor
    GroupFlowData() {}

    //! Constructor for full flow field data
    GroupFlowData(const std::string& fnbase,
                  const int          nx,
                  const int          nz,
                  const matrix       box)
    :fnbase { fnbase },
     name { "_FULL_" }
    {
        _setup_grid_and_finalize(nx, nz, box);
    }

    //! Constructor which adds `group_name` to `fnbase`
    GroupFlowData(const std::string& fnbase_original,
                  const std::string& group_name,
                  const int          nx,
                  const int          nz,
                  const matrix       box)
    :name { group_name }
    {
        fnbase.append(fnbase_original);
        fnbase.append("_");
        fnbase.append(group_name);

        _setup_grid_and_finalize(nx, nz, box);
    }

    //! Get the number of grid bins along x
    size_t nx() const { return shape[XX]; }
    //! Get the number of grid bins along z
    size_t nz() const { return shape[ZZ]; }

    //! Get the grid bin spacing along x
    real dx() const { return spacing[XX]; }
    //! Get the grid bin spacing along z
    real dz() const { return spacing[ZZ]; }

    //! Base for output file names (`fnbase_00001.dat`, ...)
    std::string fnbase;

    //! Name or identifier for group
    std::string name;

private:
    void _setup_grid_and_finalize(const int    nx,
                                  const int    nz,
                                  const matrix box)
    {
        shape = gmx::IVec{nx, 1, nz};
        spacing = gmx::RVec{
            box[XX][XX] / static_cast<real>(nx),
            box[YY][YY],
            box[ZZ][ZZ] / static_cast<real>(nz)
        };

        _finalize();
    }
};

struct FlowData {
    //! Whether or not to collect flow field data
    bool bDoFlowCollection = false;

    //! 2D grid data
    GroupFlowData flow_field;
    std::vector<GroupFlowData> group_data; // Similar data for all separate atom groups

    //! Collect flow field data at step multiples of this
    uint64_t step_collect;
    //! Average and output flow field data at step multiples of this
    uint64_t step_output;
    //! Number of samples since last output of flow field
    uint64_t num_samples;

    //! Volume of bins in grid
    double bin_volume;

    //! Empty constructor, turns off flow field collection
    FlowData() {}

    //! Constructor for full flow field
    FlowData(const std::string              &fnbase,
             const std::vector<std::string> &group_names,
             const size_t                    nx,
             const size_t                    nz,
             const matrix                    box,
             const uint64_t                  step_collect,
             const uint64_t                  step_output)
    :bDoFlowCollection { true },
     step_collect { step_collect },
     step_output { step_output },
     num_samples { 0 }
     {
        flow_field = GroupFlowData(fnbase, nx, nz, box);
        for (const auto& name : group_names)
        {
            group_data.push_back(GroupFlowData(fnbase, name, nx, nz, box));
        }

        const auto& spacing = flow_field.spacing;
        bin_volume = spacing[XX] * spacing[YY] * spacing[ZZ];
     }

    //! Zero all data for all collected flow fields and reset sample counter
    void reset_data() {
        for (auto& bin : flow_field.values)
        {
            bin.fill(0.0);
        }

        for (auto& group : group_data)
        {
            for (auto& bin : group.values)
            {
                bin.fill(0.0);
            }
        }

        num_samples = 0;
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
