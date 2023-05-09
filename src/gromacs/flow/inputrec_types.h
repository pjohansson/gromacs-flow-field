#ifndef MD_FLOW_FIELD_IR_TYPES
#define MD_FLOW_FIELD_IR_TYPES

#include <cstddef>
#include "gromacs/utility/real.h"

namespace flow
{

enum class RnemdStrategy : int {
    Edges,
    EdgeCenter,
    Count,
    Default = Edges
};
const char* enumValueToString(RnemdStrategy enumValue);

enum class RnemdAreaDefAxis : int {
    X,
    Y,
    Z,
    Count,
    Default = Z
};
const char* enumValueToString(RnemdAreaDefAxis enumValue);
size_t rnemdAxis2Index(const RnemdAreaDefAxis value);

enum class RnemdEnergyExchangeAxis : int {
    X,
    Y,
    Z,
    Count,
    Default = X
};
const char* enumValueToString(RnemdEnergyExchangeAxis enumValue);
size_t rnemdAxis2Index(const RnemdEnergyExchangeAxis value);

//! Reverse non-equilibrium molecular dynamics options
struct RNEMDOptions {
    //! Whether or not to exchange energies between areas
    bool bDoExchange = false;

    //! Axis along which exchange groups are defined
    RnemdAreaDefAxis area_def_axis;

    //! Axis for which velocity is exchanged
    RnemdEnergyExchangeAxis energy_exchange_axis;

    //! Strategy for defining exchange areas
    RnemdStrategy strategy;

    //! How often to perform the exchange (in ps)
    real tau;

    //! Size of exchange area along the defining axis
    real area_size;

    //! Additional shift of exchange areas away from the edges
    real zadj;

    //! Reference velocity, only exchange when area velocity is lower
    real ref_velocity;
};

} // namespace flow

#endif // MD_FLOW_FIELD_IR_TYPES
