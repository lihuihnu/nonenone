/**
 * @file fluid_factory.hpp
 * @brief 由强类型算例配置构造 Natural 流体系统。
 */
#pragma once

#include <common/units.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/fluid_system.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/cubic_plus_association.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>

namespace MPMC::cases
{

namespace detail
{

template <class FluidConfig, class = void>
struct ThermodynamicModelSelector
{
    static constexpr auto value = MPMC::CubicThermodynamicModel::PengRobinson;
};

template <class FluidConfig>
struct ThermodynamicModelSelector<
    FluidConfig,
    std::void_t<decltype(FluidConfig::thermodynamicModel)>>
{
    static constexpr auto value = FluidConfig::thermodynamicModel;
};

template <class FluidConfig>
inline constexpr auto thermodynamicModel =
    ThermodynamicModelSelector<FluidConfig>::value;

template <class FluidConfig, class = void>
struct HasCubicBinaryInteractionFunction : std::false_type
{
};

template <class FluidConfig>
struct HasCubicBinaryInteractionFunction<
    FluidConfig,
    std::void_t<decltype(FluidConfig::cubicBinaryInteractionCoefficient(
        0, 0, 300.0))>> : std::true_type
{
};

template <class FluidConfig, class = void>
struct HasLinearTemperatureBip : std::false_type
{
};

template <class FluidConfig>
struct HasLinearTemperatureBip<
    FluidConfig,
    std::void_t<
        decltype(FluidConfig::binaryInteractionTemperatureSlope),
        decltype(FluidConfig::binaryInteractionReferenceTemperature)>> : std::true_type
{
};

template <class FluidConfig, class = void>
struct HasVolumeTranslation : std::false_type
{
};

template <class FluidConfig>
struct HasVolumeTranslation<
    FluidConfig,
    std::void_t<decltype(FluidConfig::componentVolumeTranslation)>> : std::true_type
{
};

template <class FluidConfig, class = void>
struct HasSoreideWhitsonAqueousVolumeTranslation : std::false_type
{
};

template <class FluidConfig>
struct HasSoreideWhitsonAqueousVolumeTranslation<
    FluidConfig,
    std::void_t<decltype(FluidConfig::soreideWhitsonAqueousVolumeTranslation)>>
    : std::true_type
{
};

template <class FluidConfig, class = void>
struct HasExplicitCpaParameters : std::false_type
{
};

template <class FluidConfig>
struct HasExplicitCpaParameters<
    FluidConfig,
    std::void_t<
        decltype(FluidConfig::cpaA0),
        decltype(FluidConfig::cpaB),
        decltype(FluidConfig::cpaC1),
        decltype(FluidConfig::cpaAssociationEnergy),
        decltype(FluidConfig::cpaAssociationVolume),
        decltype(FluidConfig::cpaDonorSites),
        decltype(FluidConfig::cpaAcceptorSites)>> : std::true_type
{
};

template <class FluidConfig, class = void>
struct HasCpaCrossAssociation : std::false_type
{
};

template <class FluidConfig>
struct HasCpaCrossAssociation<
    FluidConfig,
    std::void_t<
        decltype(FluidConfig::cpaCrossAssociationEnergy),
        decltype(FluidConfig::cpaCrossAssociationVolume)>> : std::true_type
{
};

template <class FluidConfig, class = void>
struct ThermodynamicProfilerSelector
{
    static constexpr bool value = false;
};

template <class FluidConfig>
struct ThermodynamicProfilerSelector<
    FluidConfig,
    std::void_t<decltype(FluidConfig::enableThermodynamicProfiler)>>
{
    static constexpr bool value = FluidConfig::enableThermodynamicProfiler;
};

template <class FluidConfig, class = void>
struct CpaRadialDistributionSelector
{
    static constexpr auto value = MPMC::CpaRadialDistribution::Simplified;
};

template <class FluidConfig>
struct CpaRadialDistributionSelector<
    FluidConfig,
    std::void_t<decltype(FluidConfig::cpaRadialDistribution)>>
{
    static constexpr auto value = FluidConfig::cpaRadialDistribution;
};

template <class FluidConfig, class = void>
struct CpaCubicPhysicalTermSelector
{
    static constexpr auto value =
        MPMC::CpaCubicPhysicalTerm::SoaveRedlichKwong;
};

template <class FluidConfig>
struct CpaCubicPhysicalTermSelector<
    FluidConfig,
    std::void_t<decltype(FluidConfig::cpaCubicPhysicalTerm)>>
{
    static constexpr auto value = FluidConfig::cpaCubicPhysicalTerm;
};

template <class FluidConfig, class = void>
struct HasIapwsGarciaAqueousVolume : std::false_type
{
};

template <class FluidConfig>
struct HasIapwsGarciaAqueousVolume<
    FluidConfig,
    std::void_t<
        decltype(FluidConfig::useIapwsGarciaAqueousVolume),
        decltype(FluidConfig::co2Component),
        decltype(FluidConfig::aqueousVolumeWaterMolarMass),
        decltype(FluidConfig::aqueousVolumeMaximumUnsupportedMoleFraction)>>
    : std::bool_constant<FluidConfig::useIapwsGarciaAqueousVolume>
{
};

template <class FluidConfig, class = void>
struct HasMcBrideWrightAqueousViscosity : std::false_type
{
};

template <class FluidConfig>
struct HasMcBrideWrightAqueousViscosity<
    FluidConfig,
    std::void_t<
        decltype(FluidConfig::useMcBrideWrightAqueousViscosity),
        decltype(FluidConfig::co2Component),
        decltype(FluidConfig::aqueousViscosityMaximumUnsupportedMoleFraction)>>
    : std::bool_constant<FluidConfig::useMcBrideWrightAqueousViscosity>
{
};

template <class FluidConfig, class = void>
struct HasIapws2008AqueousViscosity : std::false_type
{
};

template <class FluidConfig>
struct HasIapws2008AqueousViscosity<
    FluidConfig,
    std::void_t<
        decltype(FluidConfig::useIapws2008AqueousViscosity),
        decltype(FluidConfig::aqueousViscosityMaximumSoluteMoleFraction)>>
    : std::bool_constant<FluidConfig::useIapws2008AqueousViscosity>
{
};

} // namespace detail

/**
 * @brief 从强类型算例配置构造 Natural FluidSystem。
 *
 * Factory 刻意保持 PETSc 无关，使单元测试/预检与正式算例使用完全相同的流体
 * 配置。全组分 O/G/W 模型还会登记 H2O 组分以生成富水稳定性初值；旧算例
 * 继续沿用既有水相溶解/吸附配置分支。
 */
template <class Indices, class Config>
MPMC::FluidSystem<Indices> makeFluidSystem()
{
    using Mixture = MPMC::CompositionalMixture<Indices>;
    typename MPMC::FluidSystem<Indices>::ComponentNames names{};
    for (std::size_t c = 0; c < names.size(); ++c)
        names[c] = Config::Fluid::componentNames[c];

    Mixture mixture(
        Config::Fluid::criticalTemperature,
        Config::Fluid::criticalPressure,
        Config::Fluid::criticalVolume,
        Config::Fluid::acentricFactor,
        Config::Fluid::molarMass,
        Config::Fluid::binaryInteraction);

    MPMC::CubicEquationOfState<Indices> eos(
        Config::Fluid::eosOmegaA,
        Config::Fluid::eosOmegaB,
        std::move(mixture),
        Config::Fluid::eosModelFlag,
        Config::Fluid::eosU,
        Config::Fluid::eosW);

    // 接口：旧算例继续使用常数 binaryInteraction；新算例可提供任意 k_ij(T)
    // 回调，或给出参考温度附近的线性温度斜率。
    if constexpr (detail::HasCubicBinaryInteractionFunction<typename Config::Fluid>::value)
    {
        eos.configureBinaryInteractionFunction(
            [](int i, int j, double temperature) {
                return Config::Fluid::cubicBinaryInteractionCoefficient(
                    i, j, temperature);
            });
    }
    else if constexpr (detail::HasLinearTemperatureBip<typename Config::Fluid>::value)
    {
        eos.configureBinaryInteractionFunction(
            [](int i, int j, double temperature) {
                const auto ii = static_cast<std::size_t>(i);
                const auto jj = static_cast<std::size_t>(j);
                return Config::Fluid::binaryInteraction[ii][jj]
                    + Config::Fluid::binaryInteractionTemperatureSlope[ii][jj]
                    * (temperature - Config::Fluid::binaryInteractionReferenceTemperature);
            });
    }

    if constexpr (detail::HasVolumeTranslation<typename Config::Fluid>::value)
        eos.configureVolumeTranslation(Config::Fluid::componentVolumeTranslation);

    if constexpr (detail::HasIapwsGarciaAqueousVolume<typename Config::Fluid>::value)
    {
        typename MPMC::CubicEquationOfState<Indices>::AqueousVolumeOptions aqueous;
        aqueous.waterComponent = Config::Fluid::waterComponent;
        aqueous.co2Component = Config::Fluid::co2Component;
        aqueous.waterMolarMass = Config::Fluid::aqueousVolumeWaterMolarMass;
        aqueous.maximumUnsupportedMoleFraction =
            Config::Fluid::aqueousVolumeMaximumUnsupportedMoleFraction;
        eos.configureAqueousVolume(aqueous);
    }

    if constexpr (
        detail::HasIapws2008AqueousViscosity<typename Config::Fluid>::value)
    {
        eos.configureAqueousCompositionDomain(
            Config::Fluid::waterComponent,
            Config::Fluid::aqueousViscosityMaximumSoluteMoleFraction);
    }

    if constexpr (
        detail::thermodynamicModel<typename Config::Fluid> ==
        MPMC::CubicThermodynamicModel::SoreideWhitson)
    {
        typename MPMC::CubicEquationOfState<Indices>::SoreideWhitsonOptions sw;
        sw.waterComponent = Config::Fluid::waterComponent;
        sw.salinityMolality = Config::Fluid::soreideWhitsonSalinityMolality;
        for (int component = 0; component < Indices::numComponents; ++component)
        {
            sw.aqueousWaterBip[static_cast<std::size_t>(component)] =
                [component](double temperature, double salinityMolality) {
                    return Config::Fluid::soreideWhitsonAqueousWaterBip(
                        component, temperature, salinityMolality);
                };
        }
        if constexpr (
            detail::HasSoreideWhitsonAqueousVolumeTranslation<
                typename Config::Fluid>::value)
        {
            sw.aqueousVolumeTranslation =
                Config::Fluid::soreideWhitsonAqueousVolumeTranslation;
        }
        eos.configureSoreideWhitson(std::move(sw));
    }
    else if constexpr (
        detail::thermodynamicModel<typename Config::Fluid> ==
        MPMC::CubicThermodynamicModel::CubicPlusAssociation)
    {
        using Eos = MPMC::CubicEquationOfState<Indices>;
        typename Eos::CubicPlusAssociationOptions cpa;
        constexpr double R = MPMC::units::gasConstant;

        // CPA：非缔合组分默认沿用标准 SRK 纯组分参数，使现有组分算例切换到
        // CPA 时无需在 case config 中重复填写全部 cubic 参数。
        for (int component = 0; component < Indices::numComponents; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            const double Tc = Config::Fluid::criticalTemperature[c];
            const double Pc = Config::Fluid::criticalPressure[c];
            const double omega = Config::Fluid::acentricFactor[c];
            if constexpr (
                detail::CpaCubicPhysicalTermSelector<typename Config::Fluid>::value ==
                MPMC::CpaCubicPhysicalTerm::PengRobinson)
            {
                cpa.a0[c] = 0.45724 * R * R * Tc * Tc / Pc;
                cpa.b[c] = 0.07780 * R * Tc / Pc;
                cpa.c1[c] =
                    0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
                if (Config::Fluid::eosModelFlag == 5 && omega > 0.49)
                {
                    cpa.c1[c] = 0.379642 + 1.48503 * omega
                        - 0.164423 * omega * omega
                        + 0.016666 * omega * omega * omega;
                }
            }
            else
            {
                cpa.a0[c] = 0.42748 * R * R * Tc * Tc / Pc;
                cpa.b[c] = 0.08664 * R * Tc / Pc;
                cpa.c1[c] = 0.480 + 1.574 * omega - 0.176 * omega * omega;
            }
        }

        // CPA：水默认采用常用 4C 参数化；算例仍可在下方显式覆盖全部 CPA
        // 纯组分参数。
        if constexpr (Indices::hasWater)
        {
            const std::size_t w = static_cast<std::size_t>(Config::Fluid::waterComponent);
            cpa.a0[w] = MPMC::StandardCpaWater4C::a0;
            cpa.b[w] = MPMC::StandardCpaWater4C::b;
            cpa.c1[w] = MPMC::StandardCpaWater4C::c1;
            cpa.associationEnergy[w] = MPMC::StandardCpaWater4C::epsilon;
            cpa.associationVolume[w] = MPMC::StandardCpaWater4C::beta;
            cpa.donorSites[w] = MPMC::StandardCpaWater4C::donorSites;
            cpa.acceptorSites[w] = MPMC::StandardCpaWater4C::acceptorSites;
        }

        if constexpr (detail::HasExplicitCpaParameters<typename Config::Fluid>::value)
        {
            cpa.a0 = Config::Fluid::cpaA0;
            cpa.b = Config::Fluid::cpaB;
            cpa.c1 = Config::Fluid::cpaC1;
            cpa.associationEnergy = Config::Fluid::cpaAssociationEnergy;
            cpa.associationVolume = Config::Fluid::cpaAssociationVolume;
            cpa.donorSites = Config::Fluid::cpaDonorSites;
            cpa.acceptorSites = Config::Fluid::cpaAcceptorSites;
        }
        if constexpr (detail::HasCpaCrossAssociation<typename Config::Fluid>::value)
        {
            cpa.crossAssociationEnergy = Config::Fluid::cpaCrossAssociationEnergy;
            cpa.crossAssociationVolume = Config::Fluid::cpaCrossAssociationVolume;
        }
        cpa.radialDistribution =
            detail::CpaRadialDistributionSelector<typename Config::Fluid>::value;
        cpa.physicalTerm =
            detail::CpaCubicPhysicalTermSelector<typename Config::Fluid>::value;
        eos.configureCubicPlusAssociation(std::move(cpa));
        // Natural 储层算例为等温模型。这里只缓存仅依赖温度的 CPA 纯组分/
        // 交叉参数；组成相关 mixing 与全部 AD 导数仍在每次调用时实时计算。
        eos.configureCpaTemperatureCache(Config::Fluid::temperature);
        eos.setThermodynamicProfilerEnabled(
            detail::ThermodynamicProfilerSelector<typename Config::Fluid>::value);
    }

    typename MPMC::FluidSystem<Indices>::PhaseNames phases{};
    phases[Indices::Phase::liquid] = "Oil";
    phases[Indices::Phase::vapor] = "Gas";
    if constexpr (Indices::hasWater)
        phases[Indices::Phase::water] = "Water";

    typename MPMC::FluidSystem<Indices>::PhaseArray surfaceDensity{};
    typename MPMC::FluidSystem<Indices>::PhaseArray viscosity{};
    surfaceDensity[Indices::Phase::liquid] = Config::Fluid::surfaceDensity[0];
    surfaceDensity[Indices::Phase::vapor] = Config::Fluid::surfaceDensity[1];
    viscosity[Indices::Phase::liquid] = Config::Fluid::viscosity[0];
    viscosity[Indices::Phase::vapor] = Config::Fluid::viscosity[1];
    if constexpr (Indices::hasWater)
    {
        surfaceDensity[Indices::Phase::water] = Config::Fluid::surfaceDensity[2];
        viscosity[Indices::Phase::water] = Config::Fluid::viscosity[2];
    }

    MPMC::FluidSystem<Indices> fluid(
        surfaceDensity,
        viscosity,
        std::move(eos),
        names,
        std::move(phases),
        Config::Fluid::temperature);

    using Eval = typename Indices::ValueType;
    fluid.gasRelativePermeability = [](Eval s) { return s * s; };
    fluid.waterRelativePermeability = [](Eval s) { return s * s; };
    fluid.threePhaseOilRelativePermeability = [](Eval, Eval so, Eval) { return so * so; };

    if constexpr (!Indices::fullyCompositionalThreePhase)
    {
        fluid.waterFormationVolumeFactor = [](Eval) {
            return Eval(Config::Fluid::waterFormationVolumeFactor);
        };
        fluid.waterViscosity = [](Eval) {
            return Eval(Config::Fluid::waterViscosity);
        };
    }
    else
    {
        fluid.configureFullyCompositionalThreePhase(Config::Fluid::waterComponent);
    }

    if constexpr (
        detail::HasMcBrideWrightAqueousViscosity<typename Config::Fluid>::value)
    {
        fluid.configureMcBrideWrightAqueousViscosity(
            Config::Fluid::waterComponent,
            Config::Fluid::co2Component,
            Config::Fluid::aqueousViscosityMaximumUnsupportedMoleFraction);
    }

    if constexpr (
        detail::HasIapws2008AqueousViscosity<typename Config::Fluid>::value)
    {
        fluid.configureIapws2008AqueousViscosity(
            Config::Fluid::waterComponent,
            Config::Fluid::aqueousViscosityMaximumSoluteMoleFraction);
    }

    if constexpr (Indices::hasAqueousCO2Dissolution)
    {
        fluid.configureAqueousCO2(
            Config::Dissolution::component,
            Config::Dissolution::waterMolarMass,
            Config::Dissolution::salinityMolality);
    }

    if constexpr (Indices::hasAdsorption)
    {
        fluid.configureAdsorption(
            Config::Adsorption::thetaMax,
            Config::Adsorption::coefficient);
    }

    return fluid;
}

} // namespace MPMC::cases
