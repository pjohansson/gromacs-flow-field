#ifndef MD_FLOW_FIELD
#define MD_FLOW_FIELD

#include <array>
#include <string>
#include <vector>

#include "grid.h"

struct t_filenm;
struct t_inputrec;
struct SimulationGroups;
struct t_state;
struct t_commrec;
struct t_mdatoms;

namespace gmx
{
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
enum FlowVariable : size_t
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
//!
//! We subclass `Grid3d` in order to deal with the grid bookkeeping.
//! Since the flow field data is more or less an extension of that
//! class we do this rather than creating a `Grid3d` member variable,
//! which makes working with this class more painful and harder to
//! understand.
class FlowField : public Grid3d<Bin>
{
public:
    //! Empty constructor
    FlowField() {}

    //! Constructor for full flow field data
    FlowField(const std::string& fnbase, int nx, int nz, const matrix box);

    //! Constructor which adds `group_name` to `fnbase`
    FlowField(const std::string& fnbase_original, const std::string& group_name, int nx, int nz, const matrix box);

    //! Fast method for getting the bin index for a 2D position
    //!
    //! Implemented specifically here since I have not figured out how
    //! this interface could look like in `Grid3d`. We want as fast access
    //! to the storage as possible and thus make some optimizations:
    //!
    //! 1) We know that the grid covers the entire system and not just
    //!    a subset of it. This means we can more easily calculate
    //!    the indexing and won't need to check for saturation at the
    //!    edges as the built in methods currently do.
    //!
    //! 2) After calculating the indices along each axis we use them
    //!    directly to access the bin, which skips a check for each
    //!    dimension inside the index accessor.
    //!
    //! The risk is that we make a mistake in the indexing or PBC removal
    //! which results in out-of-memory access, but this calculation is
    //! relatively easy to check.
    //!
    //! TODO: Write tests.
    size_t index_from_pos_2d(real x, real z) const;

    //! Get the number of grid bins along x
    size_t nx() const { return shape()[XX]; }
    //! Get the number of grid bins along z
    size_t nz() const { return shape()[ZZ]; }

    //! Get the grid bin spacing along x
    real dx() const { return spacing()[XX]; }
    //! Get the grid bin spacing along z
    real dz() const { return spacing()[ZZ]; }

    //! Base for output file names (`fnbase_00001.dat`, ...)
    std::string fnbase;

    //! Name or identifier for group
    std::string name;

private:
    void _setup_grid_and_finalize(int nx, int nz, const matrix box);
};


struct FlowData
{
    //! Whether or not to collect flow field data
    bool bDoFlowCollection = false;

    //! 2d flow field grid data for all groups (combined field)
    FlowField flow_field;

    //! 2d flow field data for individual groups, if multiple are selected
    std::vector<FlowField> group_data;

    //! Collect flow field data at step multiples of this
    uint64_t step_collect;

    //! Average and output flow field data at step multiples of this
    uint64_t step_output;

    //! Number of samples since last output of flow field
    uint64_t num_samples;

    //! Empty constructor, turns off flow field collection
    FlowData() {}

    //! Constructor for full flow field, turns on collection
    FlowData(const std::string&              fnbase,
             const std::vector<std::string>& group_names,
             size_t                          nx,
             size_t                          nz,
             const matrix                    box,
             uint64_t                        step_collect,
             uint64_t                        step_output);

    //! Zero all data for all collected flow fields and reset sample counter
    void reset_data();
};

//! Prepare and return a container for flow field data
FlowData init_flow_container(int                     nfile,
                             const t_filenm          fnm[],
                             const t_inputrec*       ir,
                             const SimulationGroups* groups,
                             const t_state*          state);


//! Write information about the flow field collection
void print_flow_collection_information(const FlowData& flowcr, double dt, const gmx::MDLogger& mdlog);


//! If at a collection or output step, perform actions
void flow_collect_or_output(FlowData&               flowcr,
                            int64_t                 step,
                            const t_commrec*        cr,
                            const t_inputrec*       ir,
                            const t_mdatoms*        mdatoms,
                            const t_state*          state,
                            const SimulationGroups* groups);

} // namespace flow
} // namespace gmx

#endif // MD_FLOW_FIELD
