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

/*! \brief Declarations of flow field module utilities
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

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
