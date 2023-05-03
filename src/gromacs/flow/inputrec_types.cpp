#include "gmxpre.h"
#include "inputrec_types.h"
#include "gromacs/utility/enumerationhelpers.h"

using namespace flow;

namespace flow {

/* [PETTER] Shear velocity coupling options */
const char* enumValueToString(ShearAxis_axis enumValue)
{
    static constexpr gmx::EnumerationArray<ShearAxis_axis, const char*> shear_axis_axis_names = {
        "X", "Y", "Z"
    };
    return shear_axis_axis_names[enumValue];
}

const char* enumValueToString(ShearAxis_direction enumValue)
{
    static constexpr gmx::EnumerationArray<ShearAxis_direction, const char*> shear_axis_direction_names = {
        "X", "Y", "Z"
    };
    return shear_axis_direction_names[enumValue];
}

const char* enumValueToString(ShearCouplStrategy enumValue)
{
    static constexpr gmx::EnumerationArray<ShearCouplStrategy, const char*> shear_axis_strategy_names = {
        "Edges", "Edge-Center"
    };
    return shear_axis_strategy_names[enumValue];
}

} // namespace flow
