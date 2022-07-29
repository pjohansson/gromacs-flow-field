#include "gromacs/utility/enumerationhelpers.h"

#include "inputrec_types.h"

const char* enumValueToString(flow::GridAxes enumValue)
{
    static constexpr gmx::EnumerationArray<flow::GridAxes, const char*> names = {
        "XY", "XZ", "YZ"
    };
    return names[enumValue];
}
