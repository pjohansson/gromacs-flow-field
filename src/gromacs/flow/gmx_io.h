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

/*! \brief Declares serializing routines for flow field collection.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FLOW_GMX_IO_H
#define GMX_FLOW_GMX_IO_H

#include "gmxpre.h"

#include <cstdio>

#include <vector>

struct t_inpfile;
struct t_inputrec;
struct WarningHandler;

namespace gmx
{

struct ISerializer;

namespace flow
{

struct FlowFieldOptions;

//! Read/write flow field options from an input parameter file
//
// Called during pre-processing (grompp).
//
// At the time of writing, the input file represents an .mdp file.
// The user1-grps option is overwritten with our own flow-field-grps key.
//
// \param[in/out] inputFile       File to read or write options from/into.
// \param[out]    options         Flow field options to set/get options from
// \param[out]    user1GroupsName Key of User1 group to overwrite with our own key
// \param[out]    warnings        Handler to write warnings into.
void readFlowFieldMdpOptions(std::vector<t_inpfile>* inputFile,
                             FlowFieldOptions*       options,
                             char*                   user1GroupsName,
                             WarningHandler*         warnings);

//! Verify that flow field options are set and compatible.
//
// Called during pre-processing (grompp).
//
// \param[in]  inputRec Input record with options to check.
// \param[out] warnings Handler to write errors into when detected.
void checkFlowFieldMdpOptions(const t_inputrec& inputRec, WarningHandler* warnings);

//! Do .tpx input/output for flow field options
//
// Called during pre-processing (grompp) to read/write options from/to .tpr.
//
// \param[in/out] serializer Serializer interface.
// \param[in/out] options    Flow field options to set/get options from.
void doTpxFlowFieldIo(ISerializer* serializer, FlowFieldOptions& options);

//! Print flow field options to log
//
// Called at the beginning of mdrun to print simulation options to the log.
//
// \param[out] fp      File pointer to log
// \param[in]  indent  Current indentation
// \param[in]  options Flow field options to write
void printFlowFieldOptionsToLog(FILE* fp, int indent, const FlowFieldOptions& options);

} // namespace flow
} // namespace gmx

#endif // GMX_FLOW_GMX_IO_H
