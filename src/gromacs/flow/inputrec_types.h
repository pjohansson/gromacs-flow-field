#ifndef MD_FLOW_FIELD_IR_TYPES
#define MD_FLOW_FIELD_IR_TYPES

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

} // namespace flow

#endif // MD_FLOW_FIELD_IR_TYPES
