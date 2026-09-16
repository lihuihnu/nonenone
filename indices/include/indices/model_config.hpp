/**
 * @file model_config.hpp
 * @brief Natural 模型物理开关、相行为和组分数量的编译期配置。
 */
#pragma once

namespace MPMC
{

/**
 * @brief Natural 模型使用的相行为公式类型。
 *
 * `LegacyOilGasWithIndependentWater` 保留旧模型：油/气共享组分 EOS，水相独立守恒，
 * 可选用 Henry 关系描述 CO2 溶解。
 *
 * `FullyCompositionalThreePhase` 将油富集液相、气相和水富集液相视为同一组分集的
 * 三个热力学相；H2O 只是 `NumComponents` 中的普通组分，所有组分均可在三相间分配，
 * 两组独立的逸度相等条件闭合三相平衡。
 */
enum class PhaseBehaviorModel
{
    LegacyOilGasWithIndependentWater,
    FullyCompositionalThreePhase
};

/**
 * @brief compositional natural formulation 的纯编译期模型配置。
 *
 * 最后一个模板参数选择相行为模型；默认值保持旧算例的源码接口和未知量布局不变。
 *
 * @tparam NumComponents 组分数量。FullyCompositionalThreePhase 中必须把 H2O
 *         作为普通组分计入这里。
 * @tparam HasWater 是否存在水/富水相。
 * @tparam HasWellUnknown 是否包含井底压力未知量。
 * @tparam HasAqueousCO2Dissolution 旧模型的独立水相 CO2 Henry 平衡开关。
 * @tparam HasAdsorption 是否启用吸附物理。
 * @tparam HasLandTrapping 是否启用 Land 残余气滞留物理。
 * @tparam PhaseBehavior 相行为模型。
 */
template <
    int NumComponents,
    bool HasWater,
    bool HasWellUnknown,
    bool HasAqueousCO2Dissolution = false,
    bool HasAdsorption = false,
    bool HasLandTrapping = false,
    PhaseBehaviorModel PhaseBehavior =
        PhaseBehaviorModel::LegacyOilGasWithIndependentWater>
struct CompositionalModelConfig final
{
    static_assert(
        NumComponents >= 2,
        "CompositionalModelConfig requires at least two components.");

    static constexpr PhaseBehaviorModel phaseBehavior = PhaseBehavior;
    static constexpr bool fullyCompositionalThreePhase =
        phaseBehavior == PhaseBehaviorModel::FullyCompositionalThreePhase;

    static_assert(
        !fullyCompositionalThreePhase || HasWater,
        "Fully compositional three-phase mode requires an oil/gas/water-rich phase set.");

    static_assert(
        !fullyCompositionalThreePhase || !HasAqueousCO2Dissolution,
        "Fully compositional three-phase mode already partitions CO2/H2O by EOS and must not enable the legacy aqueous-CO2 equation.");

    static_assert(
        !HasAqueousCO2Dissolution || HasWater,
        "Aqueous CO2 dissolution requires an explicit water phase.");

    static constexpr int numComponents = NumComponents;
    static constexpr bool hasWater = HasWater;
    static constexpr bool hasWellUnknown = HasWellUnknown;
    static constexpr bool hasAqueousCO2Dissolution = HasAqueousCO2Dissolution;
    static constexpr bool hasAdsorption = HasAdsorption;
    static constexpr bool hasLandTrapping = HasLandTrapping;

    /** 旧的独立水守恒方程仅在 legacy 公式中存在。 */
    static constexpr bool hasIndependentWaterConservation =
        HasWater && !fullyCompositionalThreePhase;

    static constexpr int numHydrocarbonPhases = 2;
    static constexpr int numThermodynamicPhases =
        fullyCompositionalThreePhase ? 3 : 2;
    static constexpr int numPhases =
        fullyCompositionalThreePhase
            ? 3
            : numHydrocarbonPhases + (hasWater ? 1 : 0);
};

} // namespace MPMC
