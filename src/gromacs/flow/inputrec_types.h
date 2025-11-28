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

/*! \brief Declarations of flow field module types used during t_inputrec processing.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FLOW_INPUTREC_TYPES_H
#define GMX_FLOW_INPUTREC_TYPES_H

namespace gmx
{
namespace flow
{

//! Input options for flow field collection, read from .mdp file by grompp
//!
//! In addition to these options, `user1-grps` will contain all the simulation
//! groups that the flow field data will be sampled from.
struct FlowFieldOptions
{
    //! Whether or not to collect flow field data
    bool doFlowFieldCollection = false;

    //! Interval in steps for sampling flow field data
    int nstSample = 0;

    //! Interval in steps for averaging and writing flow field data to disk
    int nstOutput = 0;

    //! Number of flow field grid bins along x
    int numBinsX = 0;

    //! Number of flow field grid bins along z
    int numBinsZ = 0;
};

} // namespace flow
} // namespace gmx

#endif // GMX_FLOW_INPUTREC_TYPES_H
