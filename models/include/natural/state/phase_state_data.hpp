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
    // 记录哪些缺失相是因为进入微量饱和度 active-set deadband 而被主动删除。
    // 该状态跨 Newton、accepted step 与 adaptive rollback 持久化，用于只对这些
    // 相施加较宽的再出现滞回，而不改变普通 missing-phase appearance 物理。
    PhasePresence phaseSuppression{std::uint8_t{0}};
    std::array<double, Indices::numComponents> vaporOilEquilibriumRatio{};
    std::array<double, Indices::numComponents> waterOilEquilibriumRatio{};
    std::array<double, 3> phaseMoleFraction{1.0, 0.0, 0.0};
    std::array<double, 3> phaseCompressibility{1.0, 1.0, 1.0};
};

} // namespace MPMC
