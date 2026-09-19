/**
 * @file phase_state_data.hpp
 * @brief O/G/W 相存在性、相分率和相组成的运行时数据。
 */
#pragma once

#include <natural/phase_state.hpp>

#include <array>

namespace MPMC
{

/**
 * @brief 存储在 phase-state Vec 中的单元二级热力学状态。
 *
 * 类型同时容纳旧油/气字段和全组分 O/G/W 字段；`PetscPhaseStateCodec` 只序列化
 * `Indices` 当前 formulation 需要的部分。这样 PETSc/grid runtime 可以共用同一
 * C++ 类型，而实际 Vec 数值布局仍由 formulation 独立定义。
 */
template <class Indices>
struct PhaseStateData
{
    // 旧油/气 formulation 字段。
    HydrocarbonPhaseState phase{HydrocarbonPhaseState::TwoPhase};
    std::array<double, Indices::numComponents> equilibriumRatio{};
    std::array<double, Indices::numComponents> overallComposition{};
    double liquidMoleFraction{0.5};
    double liquidCompressibility{1.0};
    double vaporCompressibility{1.0};

    // 全组分 O/G/W formulation 字段。
    PhasePresence phasePresence = PhasePresence::all();
    // Active-set hysteresis memory. For a missing phase this records that it
    // was removed in the trace-saturation probe band and therefore needs the
    // wider reappearance margin. If stability legitimately reintroduces that
    // phase while its equilibrium saturation is still inside the same band,
    // the bit temporarily becomes an appearance hold: the positive trace phase
    // remains active until it grows above the probe band, while S<=0 can still
    // remove it. The state persists across Newton, accepted steps and rollback.
    PhasePresence phaseSuppression{std::uint8_t{0}};
    std::array<double, Indices::numComponents> vaporOilEquilibriumRatio{};
    std::array<double, Indices::numComponents> waterOilEquilibriumRatio{};
    std::array<double, 3> phaseMoleFraction{1.0, 0.0, 0.0};
    std::array<double, 3> phaseCompressibility{1.0, 1.0, 1.0};
};

} // namespace MPMC
