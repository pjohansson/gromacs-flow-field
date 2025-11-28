#ifndef GMX_FLOW_UTILS_H
#define GMX_FLOW_UTILS_H

#include <string>
#include <vector>

//! Various utilities for this module

#include "gromacs/topology/topology.h"

namespace gmx
{
namespace flow
{
/*! \brief Get the number of groups in User1

    This is slightly complicated by how Gromacs adds a "rest" group
    to the array of names if the other groups do not add up to all
    atoms in the system. Thus, we detect if the final group is called
    exactly "rest" and if so do not count it as one of the groups. */
inline std::vector<std::string> getGroupsInUser1(const SimulationGroups& groups)
{
    std::vector<std::string> groupNames;

    for (const int globalGroupIndex : groups.groups[SimulationAtomGroupType::User1])
    {
        groupNames.emplace_back(*groups.groupNames[globalGroupIndex]);
    }

    // Remove "rest" group if it is the last value
    if (!groupNames.empty() && groupNames.back() == "rest")
    {
        groupNames.pop_back();
    }

    return groupNames;
}

} // namespace flow
} // namespace gmx

#endif // GMX_FLOW_UTILS_H
