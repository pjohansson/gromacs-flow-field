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

/*! \brief Defines serializing routines for flow field collection.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmx_io.h"

#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/txtdump.h"

#include "inputrec_types.h"

namespace gmx
{
namespace flow
{

void readFlowFieldMdpOptions(std::vector<t_inpfile>* inputFile,
                             FlowFieldOptions*       options,
                             char*                   user1GroupsName,
                             WarningHandler*         warnings)
{
    printStringNewline(inputFile, "FLOW FIELD COLLECTION");

    printStringNoNewline(inputFile, "Do flow field collection: No or Yes");
    options->doFlowFieldCollection = (getEnum<Boolean>(inputFile, "flow-field", warnings) != Boolean::No);

    printStringNoNewline(inputFile, "This selects the subset of atoms for the flow field");
    printStringNoNewline(inputFile, "collection. You can select multiple groups, in which");
    printStringNoNewline(inputFile, "case the fields for all groups combined and the fields");
    printStringNoNewline(inputFile, "for all individual groups are all collected and written");
    printStringNoNewline(inputFile, "to disk. If no groups are selected, all atoms in the");
    printStringNoNewline(inputFile, "system will be used.");
    setStringEntry(inputFile, "flow-field-grps", user1GroupsName, nullptr); // replaces user1-grps with flow-field-grps

    printStringNoNewline(inputFile, "Interval in steps between sampling flow field data");
    options->nstSample = get_eint(inputFile, "flow-nstsample", 0, warnings);
    printStringNoNewline(inputFile,
                         "Interval in steps between averaging and outputting flow field data");
    options->nstOutput = get_eint(inputFile, "flow-nstoutput", 0, warnings);

    printStringNoNewline(inputFile, "Number of flow field grid bins along x");
    options->numBinsX = get_eint(inputFile, "flow-nx", 0, warnings);
    printStringNoNewline(inputFile, "Number of flow field grid bins along z");
    options->numBinsZ = get_eint(inputFile, "flow-nz", 0, warnings);
}

void checkFlowFieldMdpOptions(const t_inputrec& inputRec, WarningHandler* warnings)
{
    const FlowFieldOptions& options = inputRec.flowFieldOptions;

    if (options.numBinsX < 1)
    {
        warnings->addError("flow-nx should be >= 1");
    }
    if (options.numBinsZ < 1)
    {
        warnings->addError("flow-nz should be >= 1");
    }
    if (options.nstSample < 1)
    {
        warnings->addError("flow-nstsample should be >= 1");
    }
    if (options.nstOutput < 1)
    {
        warnings->addError("flow-nstoutput should be >= 1");
    }

    if (options.nstOutput % options.nstSample != 0)
    {
        const std::string message = gmx::formatString(
                "flow-nstoutput (%d) should be "
                "a multiple of flow-nstsample (%d)",
                options.nstOutput,
                options.nstSample);

        warnings->addError(message);
    }

    if (inputRec.pressureCouplingOptions.epc != PressureCoupling::No)
    {
        const std::string message = gmx::formatString(
                "Pressure scaling and flow field collection were both turned "
                "on. Note that this scales flow field bins at every collection "
                "step, which may result in incorrectly sampled quantities "
                "if the box is changing size rapidly.");

        warnings->addWarning(message);
    }
}

void doTpxFlowFieldIo(ISerializer* serializer, FlowFieldOptions& options)
{
    serializer->doBool(&options.doFlowFieldCollection);
    serializer->doInt(&options.nstSample);
    serializer->doInt(&options.nstOutput);
    serializer->doInt(&options.numBinsX);
    serializer->doInt(&options.numBinsZ);
}

void printFlowFieldOptionsToLog(FILE* fp, int indent, const FlowFieldOptions& options)
{
    pr_str(fp, indent, "flow-field", booleanValueToString(options.doFlowFieldCollection));
    pr_int(fp, indent, "flow-nstsample", options.nstSample);
    pr_int(fp, indent, "flow-nstoutput", options.nstOutput);
    pr_int(fp, indent, "flow-nx", options.numBinsX);
    pr_int(fp, indent, "flow-nz", options.numBinsZ);
}

} // namespace flow
} // namespace gmx
