//! Various utilities for this module

#include "gromacs/topology/topology.h"

#ifndef MD_FLOW_FIELD_UTILS
#    define MD_FLOW_FIELD_UTILS

namespace flow
{
/*! \brief Get the number of groups in User1

    This is slightly complicated by how Gromacs adds a "rest" group
    to the array of names if the other groups do not add up to all
    atoms in the system. Thus, we detect if the final group is called
    exactly "rest" and if so do not count it as one of the groups. */
static size_t get_num_groups(const SimulationGroups* groups)
{
    size_t num_groups = 0;

    for (const auto global_group_index : groups->groups[SimulationAtomGroupType::User1])
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

} // namespace flow

#endif // MD_FLOW_FIELD_UTILS
