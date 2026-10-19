#ifndef AHF_NATIVE_CANDIDATE_TYPES_H
#define AHF_NATIVE_CANDIDATE_TYPES_H
#include <cstddef>
#include <cstdint>
namespace ahf_native {
constexpr int fields=35, stages=51;
struct RuntimeGeometry { double r, xyz[3], xud[2][3], xudd[2][6]; };
// All values here depend on patch angles/origin only, never on surface h.
struct StaticMetadata {
    double rho, sigma, ratio_first, ratio_second, sign, normalizer, origin[3];
    std::int32_t family, patch_number, is_plus, reserved;
};
}
#endif
