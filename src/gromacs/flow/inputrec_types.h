#ifndef MD_FLOW_FIELD_IR_TYPES
#define MD_FLOW_FIELD_IR_TYPES

#include "gromacs/utility/real.h"

namespace flow
{

enum class ShearCouplStrategy : int {
    Edges,
    EdgeCenter,
    Count,
    Default = Edges
};
const char* enumValueToString(ShearCouplStrategy enumValue);

enum class ShearAxis_axis : int {
    X,
    Y,
    Z,
    Count,
    Default = Z
};
const char* enumValueToString(ShearAxis_axis enumValue);

enum class ShearAxis_direction : int {
    X,
    Y,
    Z,
    Count,
    Default = X
};
const char* enumValueToString(ShearAxis_direction enumValue);

//! Reverse non-equilibrium molecular dynamics options
struct RNEMDOptions {
    //! Whether or not to exchange energies between areas
    bool bDoExchange = false;

    //! Axis along which exchange groups are defined
    ShearAxis_axis axis;

    //! Axis for which velocity is exchanged
    ShearAxis_direction direction;

    //! Strategy for defining exchange areas
    ShearCouplStrategy strategy;

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
