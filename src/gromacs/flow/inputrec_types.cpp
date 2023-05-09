#include "gmxpre.h"
#include "inputrec_types.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"

using namespace flow;

namespace flow {

/* [PETTER] Shear velocity coupling options */
const char* enumValueToString(RnemdAreaDefAxis enumValue)
{
    static constexpr gmx::EnumerationArray<RnemdAreaDefAxis, const char*> area_def_axis_names = {
        "X", "Y", "Z"
    };
    return area_def_axis_names[enumValue];
}

const char* enumValueToString(RnemdEnergyExchangeAxis enumValue)
{
    static constexpr gmx::EnumerationArray<RnemdEnergyExchangeAxis, const char*> energy_exchange_axis_names = {
        "X", "Y", "Z", "Ekin"
    };
    return energy_exchange_axis_names[enumValue];
}

const char* enumValueToString(RnemdStrategy enumValue)
{
    static constexpr gmx::EnumerationArray<RnemdStrategy, const char*> rnemd_strategy_names = {
        "Edges", "Edge-Center"
    };
    return rnemd_strategy_names[enumValue];
}

size_t rnemdAxis2Index(const RnemdAreaDefAxis value)
{
    switch (static_cast<int>(value))
    {
        case XX:
        case YY:
        case ZZ:
            return static_cast<size_t>(value);
            break;
        default:
            gmx_fatal(FARGS, "invalid value in RnemdAxis enum = %d", value);
            break;
    }
}

size_t rnemdAxis2Index(const RnemdEnergyExchangeAxis value)
{
    switch (static_cast<int>(value))
    {
        case XX:
        case YY:
        case ZZ:
            return static_cast<size_t>(value);
            break;
        default:
            gmx_fatal(FARGS, "invalid value in RnemdAxis enum = %d", value);
            break;
    }
}

} // namespace flow
