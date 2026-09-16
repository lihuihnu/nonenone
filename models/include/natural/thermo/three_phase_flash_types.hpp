/**
 * @file three_phase_flash_types.hpp
 * @brief 三相 Flash 公共 options/result 类型。
 */
#pragma once

#include <natural/thermo/phase_behavior.hpp>

#include <array>

namespace MPMC
{

/**
 * @brief 三相 flash/稳定性求解器的数值控制参数。
 *
 * 数值：当前实现组合了带保护的广义 Rachford–Rice、若干次 SSI 全局化、
 * 以及 log-K/稳定性变量上的小型稠密 Newton + 线搜索。公共结果类型不绑定具体 EOS，
 * 因此更换热力学后端时无需修改 Natural 流动方程。
 */
struct ThreePhaseFlashOptions
{
    int waterComponent{-1};
    int maximumIterations{120};
    int maximumRachfordRiceIterations{80};
    int maximumStabilityIterations{80};
    int stabilitySsiIterationsBeforeNewton{8};
    int maximumStabilityNewtonIterations{12};
    int twoPhaseSsiIterationsBeforeNewton{8};
    int maximumTwoPhaseNewtonIterations{12};
    int threePhaseSsiIterationsBeforeNewton{10};
    int maximumThreePhaseNewtonIterations{12};
    double fugacityTolerance{1.0e-9};
    double rachfordRiceTolerance{1.0e-12};
    double stabilityTolerance{1.0e-9};
    double phaseFractionTolerance{1.0e-10};
    double coincidentPhaseCompositionTolerance{1.0e-8};
    double compositionFloor{1.0e-30};
    double logKStepLimit{2.0};
    double numericalJacobianStep{1.0e-6};
};

template <class Indices>
struct ThreePhaseFlashResult
{
    static constexpr int N = Indices::numComponents;
    using Composition = std::array<double, N>;

    bool converged{false};
    int iterations{0};
    PhasePresence presence = PhasePresence::all();
    std::array<double, 3> phaseMoleFraction{1.0, 0.0, 0.0};
    std::array<double, 3> saturation{1.0, 0.0, 0.0};
    std::array<Composition, 3> composition{};
    std::array<double, 3> compressibility{1.0, 1.0, 1.0};
    std::array<double, 3> molarDensity{0.0, 0.0, 0.0};
    Composition vaporOilK{};
    Composition waterOilK{};
};

template <class Indices>
struct ThreePhaseStabilityResult
{
    static constexpr int N = Indices::numComponents;
    using Composition = std::array<double, N>;

    bool stable{true};
    bool valid{true};
    PhasePresence testedPresence{};
    std::array<bool, 3> missingPhaseUnstable{false, false, false};
    std::array<Composition, 3> incipientComposition{};
    std::array<double, 3> trialSum{1.0, 1.0, 1.0};
};


} // namespace MPMC
