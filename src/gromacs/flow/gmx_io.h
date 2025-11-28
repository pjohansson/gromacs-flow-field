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
