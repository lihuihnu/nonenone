/**
 * @file numerics.hpp
 * @brief Natural 热力学与非线性求解共用的数值阈值。
 */
#pragma once

namespace MPMC
{

/**
 * @brief Natural formulation 的数值常量。
 *
 * 全部 Natural 核心共用这一组常量，避免同一物理模型出现多套隐式数值路径。
 */
struct NaturalNumerics final
{
    static constexpr double minimumComposition = 1.0e-8;
    static constexpr double minimumSaturation = 1.0e-8;
    static constexpr double phaseAppearanceSaturation = 1.0e-6;

    /**
     * @brief 全组分三相 active-set 的热力学相边界探测带。
     *
     * 该值不是“强制删相”的物理阈值，更不是残余饱和度。活动相进入该带后，
     * 程序只获得一次构造约化相集的资格，随后仍必须通过 restricted flash 和
     * missing-phase stability 检查。把探测带设得略宽于机器级 trace saturation，
     * 可以避免 Newton 在相边界外侧因退化的 disappearing-phase 变量而渐近停滞；
     * 热力学上明确稳定存在的相会立即由 stability 检查恢复。
     */
    static constexpr double phaseBoundaryProbeSaturation = 1.0e-4;

    // Backward-compatible name retained for downstream case code.  New code
    // should use phaseBoundaryProbeSaturation to make clear that this is only
    // a thermodynamic active-set probe, not a physical residual saturation.
    static constexpr double phaseDisappearanceSaturation = phaseBoundaryProbeSaturation;

    /**
     * @brief 普通缺失相重新出现所需的稳定性超额量。
     *
     * 对初始就缺失、或并非由数值相边界主动删除的相，继续保持严格的
     * stability 判据，避免把真实的新相生成延后。
     */
    static constexpr double phaseAppearanceStabilityMargin = 1.0e-6;

    /**
     * @brief 因微量饱和度主动删除的相所使用的再出现滞回裕量。
     *
     * 只有 phase-state 中记录为 suppressed-by-hysteresis 的相使用该值。它与
     * phaseBoundaryProbeSaturation 采用同一 1e-4 数量级 deadband：这样一个刚被
     * active-set 删除的 trace phase 不会因为 TPD/trial-sum 的 1e-5 级数值噪声在
     * 下一 Newton 立即重生；missing-phase stability 一旦明确越过该裕量，相仍会
     * 正常重新出现。
     */
    static constexpr double phaseHysteresisReappearanceMargin = 1.0e-4;
    static constexpr double minimumNormalizationDenominator = 1.0e-14;
    static constexpr double phaseComparisonTolerance = 1.0e-12;

    /**
     * @brief 判定组分在某相中数值缺失的 active-set 阈值。
     *
     * The N-1 composition parameterization cannot represent a dependent component
     * far below machine precision relative to unity.  When a flash places a
     * component below this threshold in the non-reference phase, the fully
     * compositional residual uses the composition boundary condition x_i=0
     * instead of an impossible fugacity-equality row.
     */
    static constexpr double phaseEquilibriumTraceComposition = 1.0e-14;

    static constexpr double maximumSaturationNewtonChange = 0.1;
    static constexpr double maximumCompositionNewtonChange = 0.1;
    static constexpr double maximumAqueousCO2NewtonChange = 0.02;
    static constexpr double maximumRelativePressureNewtonChange = 0.25;

    static constexpr int phaseStabilityMaximumIterations = 20000;
    static constexpr double phaseStabilityTolerance = 1.0e-10;

};

} // namespace MPMC
