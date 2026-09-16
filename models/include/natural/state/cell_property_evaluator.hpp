/**
 * @file cell_property_evaluator.hpp
 * @brief 由主变量和相状态计算单元热力学与流动属性。
 */
#pragma once

#include <stdexcept>
#include <natural/fluid_system.hpp>
#include <natural/numerics.hpp>
#include <natural/physics/land_trapping.hpp>
#include <natural/state/cell_state.hpp>

#include <array>
#include <type_traits>

namespace MPMC
{

/**
 * @brief 纯单元物性计算器。
 *
 * 输入一个已经构造好的 CellState 和岩石局部量，输出完整 CellProperties。
 * 不读 Vec、不知道 cell id、不知道 StructuredGrid/CpGrid。
 */
template <class Indices>
class CellPropertyEvaluator final
{
public:
    using Scalar = typename Indices::ValueType;
    using State = CellState<Indices, Scalar>;
    using Properties = CellProperties<Indices, Scalar>;
    using Composition = std::array<Scalar, Indices::numComponents>;

    explicit CellPropertyEvaluator(const FluidSystem<Indices> &fluid)
        : fluid_(fluid)
    {
        validateConstitutiveModels_();
    }

    [[nodiscard]] Properties evaluate(
        const State &state,
        const Scalar &basePorosity,
        [[maybe_unused]] double maximumGasSaturation = 0.0,
        [[maybe_unused]] double landConstant = 2.0,
        [[maybe_unused]] const Composition *overallComposition = nullptr) const
    {
        if constexpr (Indices::fullyCompositionalThreePhase)
            return evaluateFullyCompositional_(
                state, basePorosity, maximumGasSaturation, landConstant);
        else
            return evaluateLegacy_(
                state, basePorosity, maximumGasSaturation, landConstant, overallComposition);
    }

private:
    [[nodiscard]] Properties evaluateLegacy_(
        const State &state,
        const Scalar &basePorosity,
        [[maybe_unused]] double maximumGasSaturation,
        [[maybe_unused]] double landConstant,
        const Composition *overallComposition) const
    {
        Properties properties;
        properties.porosity = basePorosity;

        const int liquid = Indices::Phase::liquid;
        const int vapor = Indices::Phase::vapor;
        properties.saturation[static_cast<std::size_t>(liquid)] = state.liquidSaturation;
        properties.saturation[static_cast<std::size_t>(vapor)] = state.vaporSaturation;

        Composition x = state.liquidMoleFraction;
        Composition y = state.vaporMoleFraction;

        if (state.hydrocarbonPhaseState == HydrocarbonPhaseState::VaporOnly)
        {
            y = x;
        }
        else if (state.hydrocarbonPhaseState == HydrocarbonPhaseState::LiquidOnly)
        {
            y = x;
            clearCompositionDerivatives(y);
        }

        Scalar liquidZ = 1.0;
        Scalar vaporZ = 1.0;
        typename CubicEquationOfState<Indices>::template PhaseResult<Scalar> liquidThermo;
        typename CubicEquationOfState<Indices>::template PhaseResult<Scalar> vaporThermo;

        if (!fluid_.eos.usesEquilibriumConstants())
        {
            if (fluid_.eos.usesCubicPlusAssociation())
            {
                liquidThermo = fluid_.eos.phaseResult(
                    state.pressure, Scalar(fluid_.temperature), x, true);
                vaporThermo = fluid_.eos.phaseResult(
                    state.pressure, Scalar(fluid_.temperature), y, false);
                liquidZ = liquidThermo.compressibility;
                vaporZ = vaporThermo.compressibility;
            }
            else
            {
                const auto params = fluid_.eos.mixingParameters(
                    state.pressure, Scalar(fluid_.temperature));
                const auto mixL = fluid_.eos.phaseMixing(x, params);
                const auto mixV = fluid_.eos.phaseMixing(y, params);
                liquidZ = fluid_.eos.liquidRoot(mixL.A, mixL.B);
                vaporZ = fluid_.eos.vaporRoot(mixV.A, mixV.B);

                applyVaporZLinearization(vaporZ);

                liquidThermo = fluid_.eos.fugacityAtRoot(
                    state.pressure, x, params, mixL, liquidZ);
                vaporThermo = fluid_.eos.fugacityAtRoot(
                    state.pressure, y, params, mixV, vaporZ);
            }

            if (state.hydrocarbonPhaseState == HydrocarbonPhaseState::LiquidOnly)
            {
                vaporZ = liquidZ;
                vaporThermo = liquidThermo;
            }
            else if (state.hydrocarbonPhaseState == HydrocarbonPhaseState::VaporOnly)
            {
                liquidZ = vaporZ;
                liquidThermo = vaporThermo;
            }

            auto [xMass, rhoL, muL] = fluid_.eosFlowProperties(
                state.pressure, x, liquidZ, Scalar(fluid_.temperature));
            auto [yMass, rhoV, muV] = fluid_.eosFlowProperties(
                state.pressure, y, vaporZ, Scalar(fluid_.temperature));

            properties.density[static_cast<std::size_t>(liquid)] = rhoL;
            properties.density[static_cast<std::size_t>(vapor)] = rhoV;
            properties.viscosity[static_cast<std::size_t>(liquid)] = muL;
            properties.viscosity[static_cast<std::size_t>(vapor)] = muV;
            properties.massFraction[static_cast<std::size_t>(liquid)] = xMass;
            properties.massFraction[static_cast<std::size_t>(vapor)] = yMass;
            properties.fugacity[0] = liquidThermo.fugacity;
            properties.fugacity[1] = vaporThermo.fugacity;
        }
        else
        {
            // 黑油/K-value 模式：液相伪逸度 fL=K*x，气相伪逸度 fV=y。
            auto [xMass, rhoL, muL] = fluid_.blackOilFlowProperties(
                state.pressure, x, true);
            auto [yMass, rhoV, muV] = fluid_.blackOilFlowProperties(
                state.pressure, y, false);
            properties.density[static_cast<std::size_t>(liquid)] = rhoL;
            properties.density[static_cast<std::size_t>(vapor)] = rhoV;
            properties.viscosity[static_cast<std::size_t>(liquid)] = muL;
            properties.viscosity[static_cast<std::size_t>(vapor)] = muV;
            properties.massFraction[static_cast<std::size_t>(liquid)] = xMass;
            properties.massFraction[static_cast<std::size_t>(vapor)] = yMass;

            // K 值可依赖总体组成 z。纯单元计算器不读取分布式 phase-state，
            // 因此由调用者显式传入 z；未提供时使用液相组成 x。
            const Composition &z = overallComposition != nullptr
                ? *overallComposition
                : x;
            const auto K = fluid_.eos.evaluateEquilibriumConstants(
                state.pressure, Scalar(fluid_.temperature), z);
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                properties.fugacity[0][c] = K[c] * x[c];
                properties.fugacity[1][c] = y[c];
            }
        }

        Scalar gasSaturationForRelperm = state.vaporSaturation;
        if (state.hydrocarbonPhaseState == HydrocarbonPhaseState::LiquidOnly)
            clearDerivatives(gasSaturationForRelperm);
        Scalar liquidSaturationForRelperm = state.liquidSaturation;
        if (state.hydrocarbonPhaseState == HydrocarbonPhaseState::VaporOnly)
            clearDerivatives(liquidSaturationForRelperm);

        if constexpr (Indices::hasLandTrapping)
        {
            const LandTrappingModel land(landConstant);
            properties.trappedGasSaturation = land.trappedGasSaturation(
                gasSaturationForRelperm, maximumGasSaturation);
        }
        properties.freeGasSaturation =
            gasSaturationForRelperm - properties.trappedGasSaturation;

        const Scalar krg = fluid_.gasRelativePermeability(properties.freeGasSaturation);
        Scalar kro = 0.0;
        if constexpr (Indices::numPhases == 3)
            kro = fluid_.threePhaseOilRelativePermeability(
                state.waterSaturation,
                liquidSaturationForRelperm,
                gasSaturationForRelperm);
        else
            kro = fluid_.oilRelativePermeability(liquidSaturationForRelperm);

        properties.mobility[static_cast<std::size_t>(liquid)] =
            kro / properties.viscosity[static_cast<std::size_t>(liquid)];
        properties.mobility[static_cast<std::size_t>(vapor)] =
            krg / properties.viscosity[static_cast<std::size_t>(vapor)];

        if constexpr (Indices::hasWater)
        {
            const int water = Indices::Phase::water;
            const Scalar krw = fluid_.waterRelativePermeability(state.waterSaturation);
            const Scalar muw = fluid_.waterViscosity(state.pressure);
            const Scalar bw = fluid_.waterFormationVolumeFactor(state.pressure);
            const Scalar rhoW = bw * fluid_.surfaceDensity[static_cast<std::size_t>(water)];

            properties.saturation[static_cast<std::size_t>(water)] = state.waterSaturation;
            properties.density[static_cast<std::size_t>(water)] = rhoW;
            properties.viscosity[static_cast<std::size_t>(water)] = muw;
            properties.mobility[static_cast<std::size_t>(water)] = krw / muw;

            if constexpr (Indices::hasAqueousCO2Dissolution)
            {
                properties.aqueousCO2MassFraction =
                    fluid_.aqueousCO2MassFraction(state.aqueousCO2MoleFraction);
                properties.aqueousCO2Fugacity =
                    fluid_.aqueousCO2Fugacity(state.pressure, state.aqueousCO2MoleFraction);
            }
        }

        if constexpr (Indices::hasAdsorption)
        {
            properties.adsorbedVolume = fluid_.adsorbedVolume(
                state.pressure, y);
        }

        if (fluid_.poreVolumeMultiplier)
            properties.porosity = fluid_.poreVolumeMultiplier(state.pressure) * properties.porosity;

        return properties;
    }

    [[nodiscard]] Properties evaluateFullyCompositional_(
        const State &state,
        const Scalar &basePorosity,
        [[maybe_unused]] double maximumGasSaturation,
        [[maybe_unused]] double landConstant) const
    {
        static_assert(Indices::numPhases == 3,
                      "Fully compositional mode requires exactly three phases.");
        if (fluid_.eos.usesEquilibriumConstants())
            throw std::logic_error(
                "Fully compositional three-phase flow requires EOS fugacities.");

        Properties properties;
        properties.porosity = basePorosity;

        std::array<Composition, 3> composition{
            state.liquidMoleFraction,
            state.vaporMoleFraction,
            state.aqueousMoleFraction};
        const std::array<Scalar, 3> saturation{
            state.liquidSaturation,
            state.vaporSaturation,
            state.waterSaturation};
        const std::array<CompositionalPhase, 3> phases{
            CompositionalPhase::Oil,
            CompositionalPhase::Gas,
            CompositionalPhase::Water};

        // Reduced phase sets retain finite placeholder compositions for absent
        // phases.  A hydrocarbon-rich overall composition is a useful generic
        // flash placeholder, but it lies outside the deliberately strict
        // IAPWS-Garcia/McBride-Wright aqueous domain.  Property evaluation still
        // visits every phase to keep fixed-size residual storage, so evaluate an
        // absent water phase at a pure-H2O reference state.  Its saturation and
        // mobility remain exactly zero; an active water phase continues to use
        // the physical flash composition and therefore remains domain checked.
        if (!state.phasePresence.contains(CompositionalPhase::Water))
        {
            const int waterComponent = fluid_.fullyCompositionalWaterComponent;
            if (waterComponent < 0 || waterComponent >= Indices::numComponents)
                throw std::logic_error(
                    "Fully compositional water phase requires a valid H2O component index.");
            Composition reference{};
            reference[static_cast<std::size_t>(waterComponent)] = Scalar(1.0);
            composition[static_cast<std::size_t>(Indices::Phase::water)] = reference;
        }

        const Scalar temperature = Scalar(fluid_.temperature);

        if (!fluid_.eos.usesCubicPlusAssociation())
        {
            // PR parameters depend on (p,T) but not phase identity.  SW uses the
            // same non-aqueous parameters for oil/gas and a separate aqueous BIP
            // set only for the water-rich phase.  Build each distinct parameter
            // set once per cell, then keep the original phaseMixing/root/fugacity
            // operations unchanged for every phase.
            const auto nonAqueousParameters = fluid_.eos.mixingParameters(
                state.pressure, temperature, CompositionalPhase::Oil);

            for (int phase = 0; phase < 3; ++phase)
            {
                const std::size_t p = static_cast<std::size_t>(phase);
                properties.saturation[p] = saturation[p];

                if (phases[p] == CompositionalPhase::Water &&
                    fluid_.eos.usesSoreideWhitson())
                {
                    const auto aqueousParameters = fluid_.eos.mixingParameters(
                        state.pressure, temperature, CompositionalPhase::Water);
                    const auto mix = fluid_.eos.phaseMixing(
                        composition[p], aqueousParameters);
                    const Scalar z = fluid_.eos.liquidRoot(mix.A, mix.B);
                    const auto thermo = fluid_.eos.fugacityAtRoot(
                        state.pressure, composition[p], aqueousParameters, mix, z);
                    const auto [mass, density, viscosity] = fluid_.eosFlowProperties(
                        state.pressure, composition[p], z, temperature, phases[p]);
                    properties.density[p] = density;
                    properties.viscosity[p] = viscosity;
                    properties.massFraction[p] = mass;
                    properties.fugacity[p] = thermo.fugacity;
                }
                else
                {
                    const auto mix = fluid_.eos.phaseMixing(
                        composition[p], nonAqueousParameters);
                    const Scalar z = phases[p] == CompositionalPhase::Gas
                        ? fluid_.eos.vaporRoot(mix.A, mix.B)
                        : fluid_.eos.liquidRoot(mix.A, mix.B);
                    const auto thermo = fluid_.eos.fugacityAtRoot(
                        state.pressure, composition[p], nonAqueousParameters, mix, z);
                    const auto [mass, density, viscosity] = fluid_.eosFlowProperties(
                        state.pressure, composition[p], z, temperature, phases[p]);
                    properties.density[p] = density;
                    properties.viscosity[p] = viscosity;
                    properties.massFraction[p] = mass;
                    properties.fugacity[p] = thermo.fugacity;
                }
            }
        }
        else
        {
            // CPA density/root evaluation is composition dependent; keep the full
            // phaseResult path for each phase.
            for (int phase = 0; phase < 3; ++phase)
            {
                const std::size_t p = static_cast<std::size_t>(phase);
                properties.saturation[p] = saturation[p];
                const auto thermo = fluid_.eos.phaseResult(
                    state.pressure, temperature, composition[p], phases[p]);
                const auto [mass, density, viscosity] = fluid_.eosFlowProperties(
                    state.pressure, composition[p], thermo.compressibility,
                    temperature, phases[p]);
                properties.density[p] = density;
                properties.viscosity[p] = viscosity;
                properties.massFraction[p] = mass;
                properties.fugacity[p] = thermo.fugacity;
            }
        }

        Scalar gasSaturationForRelperm = state.vaporSaturation;
        if (!state.phasePresence.contains(CompositionalPhase::Gas))
            clearDerivatives(gasSaturationForRelperm);

        if constexpr (Indices::hasLandTrapping)
        {
            const LandTrappingModel land(landConstant);
            properties.trappedGasSaturation = land.trappedGasSaturation(
                gasSaturationForRelperm, maximumGasSaturation);
        }
        properties.freeGasSaturation =
            gasSaturationForRelperm - properties.trappedGasSaturation;

        const Scalar kro = fluid_.threePhaseOilRelativePermeability(
            state.waterSaturation, state.liquidSaturation, gasSaturationForRelperm);
        const Scalar krg = fluid_.gasRelativePermeability(properties.freeGasSaturation);
        const Scalar krw = fluid_.waterRelativePermeability(state.waterSaturation);

        const std::array<Scalar, 3> kr{kro, krg, krw};
        for (int phase = 0; phase < 3; ++phase)
        {
            const std::size_t p = static_cast<std::size_t>(phase);
            if (!state.phasePresence.contains(phases[p]))
            {
                properties.mobility[p] = Scalar(0.0);
                continue;
            }
            properties.mobility[p] = kr[p] / properties.viscosity[p];
        }

        if constexpr (Indices::hasAdsorption)
            properties.adsorbedVolume = fluid_.adsorbedVolume(
                state.pressure, state.vaporMoleFraction);

        if (fluid_.poreVolumeMultiplier)
            properties.porosity =
                fluid_.poreVolumeMultiplier(state.pressure) * properties.porosity;
        return properties;
    }

    /**
     * @brief 在 evaluator 构造时一次性验证岩石-流体本构闭包。
     *
     * 这些 std::function 在一个模拟过程中保持不变；把配置检查移出
     * 每单元/每 Newton 次的热路径，可避免重复分支，同时让缺失配置更早失败。
     */
    void validateConstitutiveModels_() const
    {
        if (!fluid_.gasRelativePermeability)
            throw std::logic_error(
                "FluidSystem::gasRelativePermeability is not configured.");

        if constexpr (Indices::numPhases == 3)
        {
            if (!fluid_.threePhaseOilRelativePermeability)
                throw std::logic_error(
                    "FluidSystem::threePhaseOilRelativePermeability is not configured.");
        }
        else if (!fluid_.oilRelativePermeability)
        {
            throw std::logic_error(
                "FluidSystem::oilRelativePermeability is not configured.");
        }

        if constexpr (Indices::hasWater)
        {
            if (!fluid_.waterRelativePermeability)
                throw std::logic_error(
                    "FluidSystem::waterRelativePermeability is not configured.");

            if constexpr (!Indices::fullyCompositionalThreePhase)
            {
                if (!fluid_.waterViscosity || !fluid_.waterFormationVolumeFactor)
                    throw std::logic_error(
                        "Water viscosity/Bw callbacks are not fully configured.");
            }
        }
    }

    static void clearDerivatives(Scalar &value)
    {
        if constexpr (!std::is_same_v<Scalar, double>)
            value.clearDerivatives();
    }

    static void clearCompositionDerivatives(Composition &composition)
    {
        if constexpr (!std::is_same_v<Scalar, double>)
            for (auto &entry : composition)
                entry.clearDerivatives();
    }

    /**
     * @brief 固定气相压缩因子前 N 个 AD 导数槽。
     *
     * Natural EOS 线性化在求得气相 Z 后将导数槽 `[0, numComponents)` 清零，
     * 再用该 Z 计算气相逸度。这里显式封装该规则，
     * 避免在 EOS 根求解代码中直接访问 AD 内部存储。
     */
    static void applyVaporZLinearization(Scalar &z)
    {
        if constexpr (!std::is_same_v<Scalar, double>)
        {
            const int n = std::min(Indices::numComponents, z.size());
            for (int i = 0; i < n; ++i)
                z.setDerivative(i, 0.0);
        }
    }

    const FluidSystem<Indices> &fluid_;
};

} // namespace MPMC
