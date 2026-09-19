/**
 * @file fluid_system.hpp
 * @brief Natural 模型共享的流体物性、热力学和附加物理配置。
 */
#pragma once

#include <natural/physics/aqueous_co2.hpp>
#include <natural/physics/competitive_adsorption.hpp>
#include <natural/properties/black_oil_properties.hpp>
#include <natural/properties/aqueous_viscosity.hpp>
#include <natural/properties/compositional_properties.hpp>
#include <natural/thermo/cubic_eos.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace MPMC
{

/**
 * @brief Natural 模型跨网格后端共享的储层流体描述。
 *
 * 该类拥有 EOS、物性模型和算例提供的本构关系，但不依赖网格、MPI 或 PETSc；
 * 因而同一个 `FluidSystem` 可以直接用于 StructuredGrid 和 CpGrid 后端。
 *
 * 相对渗透率为无量纲；黏度单位为 Pa·s；水体积系数和孔隙体积乘子均为无量纲。
 */
template <class Indices>
class FluidSystem
{
public:
    static constexpr int numComponents = Indices::numComponents;
    static constexpr int numPhases = Indices::numPhases;

    using ValueType = typename Indices::ValueType;
    using Composition = std::array<ValueType, numComponents>;
    using PhaseArray = std::array<double, numPhases>;
    using ComponentNames = std::array<std::string, numComponents>;
    using PhaseNames = std::array<std::string, numPhases>;
    using UnaryProperty = std::function<ValueType(ValueType)>;
    using ThreePhaseRelativePermeability =
        std::function<ValueType(ValueType, ValueType, ValueType)>;

    FluidSystem() = default;

    /** @brief 由定长数组构造流体系统；适合正式算例配置。 */
    FluidSystem(
        PhaseArray surfaceDensityValues,
        PhaseArray viscosityValues,
        CubicEquationOfState<Indices> eosModel,
        ComponentNames names,
        PhaseNames phases,
        double fluidTemperature)
        : surfaceDensity(std::move(surfaceDensityValues)),
          phaseViscosity(std::move(viscosityValues)),
          eos(std::move(eosModel)),
          compositionalProperties(eos.mixture()),
          componentNames(std::move(names)),
          phaseNames(std::move(phases)),
          temperature(fluidTemperature)
    {
        validate();
    }

    /** @brief 由 initializer_list 构造流体系统；适合简短测试/示例。 */
    FluidSystem(
        std::initializer_list<double> surfaceDensityValues,
        std::initializer_list<double> viscosityValues,
        CubicEquationOfState<Indices> eosModel,
        std::initializer_list<std::string> names,
        std::initializer_list<std::string> phases,
        double fluidTemperature)
        : eos(std::move(eosModel)),
          compositionalProperties(eos.mixture()),
          temperature(fluidTemperature)
    {
        copyList(surfaceDensityValues, surfaceDensity, "surface density");
        copyList(viscosityValues, phaseViscosity, "phase viscosity");
        copyList(names, componentNames, "component name");
        copyList(phases, phaseNames, "phase name");
        validate();
    }

    /**
     * @brief 由 EOS 相组成和压缩因子计算质量分数、密度和黏度。
     *
     * `compressibility` 必须与当前相组成来自同一热力学状态，避免物性与 flash 状态不一致。
     */
    template <class Scalar>
    [[nodiscard]] std::tuple<std::array<Scalar, numComponents>, Scalar, Scalar>
    eosFlowProperties(
        const Scalar &pressure,
        const std::array<Scalar, numComponents> &moleFraction,
        const Scalar &compressibility,
        const Scalar &fluidTemperature,
        CompositionalPhase phase = CompositionalPhase::Oil) const
    {
        const auto mass = eos.massFractions(moleFraction);

        // Density and LBC viscosity need the same phase molar density.  Build it
        // once per phase instead of evaluating p/(RTZ) twice on the dominant
        // non-volume-translated path.  The EOS path remains authoritative when
        // a volume translation is configured.
        const auto rhoMolar = eos.usesVolumeTranslation()
            ? eos.molarDensity(
                  pressure, fluidTemperature, moleFraction, compressibility, phase)
            : compositionalProperties.molarDensity(
                  pressure, compressibility, fluidTemperature);
        const auto density = compositionalProperties.densityFromMolarDensity(
            moleFraction, rhoMolar);
        const auto cubicRhoMolar = compositionalProperties.molarDensity(
            pressure, compressibility, fluidTemperature);
        Scalar viscosity = compositionalProperties.viscosityFromMolarDensity(
            moleFraction, cubicRhoMolar, fluidTemperature);
        if (phase == CompositionalPhase::Water)
        {
            if (iapws2008AqueousViscosity_)
                viscosity = iapws2008AqueousViscosity_->viscosity(
                    pressure, fluidTemperature, moleFraction);
            else if (mcBrideWrightAqueousViscosity_)
                viscosity = mcBrideWrightAqueousViscosity_->viscosity(
                    pressure, fluidTemperature, moleFraction);
        }
        if (flowViscosityOverride)
        {
            Composition values{};
            for (int i=0; i<numComponents; ++i)
                values[i] = ValueType(mass[i]);
            const auto v = flowViscosityOverride(ValueType(pressure), values);
            if constexpr (std::is_same_v<Scalar, double>) viscosity = scalarValue(v);
            else viscosity = Scalar(v);
        }
        return {mass, density, viscosity};
    }

    /** @brief 对 H2O 主导的 Water 角色相启用 IAPWS-2008 黏度。 */
    void configureIapws2008AqueousViscosity(
        int waterComponent,
        double maximumSoluteMoleFraction = 0.02)
    {
        iapws2008AqueousViscosity_.emplace(
            waterComponent, maximumSoluteMoleFraction);
    }

    /** @brief 仅对 Water 相角色启用 McBride-Wright H2O-CO2 黏度闭包。 */
    void configureMcBrideWrightAqueousViscosity(
        int waterComponent,
        int co2Component,
        double maximumUnsupportedMoleFraction = 1.0e-4)
    {
        mcBrideWrightAqueousViscosity_.emplace(
            waterComponent, co2Component,
            maximumUnsupportedMoleFraction);
    }

    /** @brief 计算旧 K-value/黑油路径中单相的密度和黏度。 */
    template <class Scalar>
    [[nodiscard]] std::tuple<std::array<Scalar, numComponents>, Scalar, Scalar>
    blackOilFlowProperties(
        const Scalar &pressure,
        const std::array<Scalar, numComponents> &moleFraction,
        bool liquid) const
    {
        const auto mass = eos.massFractions(moleFraction);
        const auto density = blackOilProperties.computeMolarDensity(
            pressure, moleFraction, liquid);
        const auto viscosity = blackOilProperties.computeViscosity(
            pressure, moleFraction, liquid);
        return {mass, density, viscosity};
    }


    /**
     * @brief 指定全组分三相 flash 中用于构造富水初值的 H2O 组分。
     *
     * H2O 仍是普通 EOS 组分；该索引只帮助稳定性/flash 构造富水初值，
     * 不施加不互溶条件，也不会额外引入 Henry 方程。
     */
    void configureFullyCompositionalThreePhase(int waterComponentIndex)
    {
        if constexpr (!Indices::fullyCompositionalThreePhase)
            throw std::logic_error("Fully compositional three-phase configuration requires the matching model config.");
        if (waterComponentIndex < 0 || waterComponentIndex >= numComponents)
            throw std::invalid_argument("Invalid H2O component index.");
        fullyCompositionalWaterComponent = waterComponentIndex;
    }

    /**
     * @brief 为旧独立水相模型启用指定组分的 CO2 溶解关系。
     *
     * 该模型负责水相 CO2 摩尔分数/质量分数换算，并按温度和盐度计算溶解态逸度。
     */
    void configureAqueousCO2(
        int componentIndex,
        double waterMolarMass,
        double salinityMolality = 0.0)
    {
        if (componentIndex < 0 || componentIndex >= numComponents)
            throw std::invalid_argument("Invalid CO2 component index.");
        aqueousCO2Component = componentIndex;
        aqueousCO2_.emplace(
            temperature,
            eos.mixture().molecularWeight(componentIndex),
            waterMolarMass,
            salinityMolality);
    }

    /**
     * @brief 对全部气相组分启用竞争 Langmuir 吸附。
     *
     * 数学：`V_i = thetaMax_i b_i p y_i / (1 + sum_j b_j p y_j)`。
     */
    void configureAdsorption(
        const std::array<double, numComponents> &thetaMax,
        const std::array<double, numComponents> &coefficient)
    {
        adsorption_.emplace(thetaMax, coefficient);
    }

    /** @brief 将水相 CO2 摩尔分数换算为蓄积/输出使用的质量分数。 */
    template <class Scalar>
    [[nodiscard]] Scalar aqueousCO2MassFraction(const Scalar &x) const
    {
        return aqueousCO2Model().massFraction(x);
    }

    /** @brief 计算溶解态 CO2 逸度 `f_w = x_CO2 H(T,p,m_s)`。 */
    template <class Scalar>
    [[nodiscard]] Scalar aqueousCO2Fugacity(const Scalar &pressure, const Scalar &x) const
    {
        return aqueousCO2Model().fugacity(pressure, x);
    }

    /** @brief 返回各气相组分的竞争 Langmuir 吸附体积。 */
    template <class Scalar>
    [[nodiscard]] std::array<Scalar, numComponents> adsorbedVolume(
        const Scalar &gasPressure,
        const std::array<Scalar, numComponents> &gasMolarFractions) const
    {
        return adsorptionModel().adsorbedVolume(gasPressure, gasMolarFractions);
    }

    /** @brief 返回已配置的水相 CO2 模型；未启用时抛出异常。 */
    [[nodiscard]] const AqueousCO2Model &aqueousCO2Model() const
    {
        if (!aqueousCO2_)
            throw std::logic_error("Aqueous CO2 model has not been configured.");
        return *aqueousCO2_;
    }

    /** @brief 返回已配置的竞争吸附模型；未启用时抛出异常。 */
    [[nodiscard]] const CompetitiveLangmuirAdsorption<numComponents> &adsorptionModel() const
    {
        if (!adsorption_)
            throw std::logic_error("Adsorption model has not been configured.");
        return *adsorption_;
    }

    PhaseArray surfaceDensity{};       ///< Surface/reference density, kg/m3.
    PhaseArray phaseViscosity{};       ///< Constant fallback phase viscosity, Pa.s.
    CubicEquationOfState<Indices> eos{}; ///< EOS, flash and fugacity calculations.
    CompositionalPropertyModel<Indices> compositionalProperties{}; ///< EOS density/viscosity closure.
    BlackOilPropertyModel<Indices> blackOilProperties{}; ///< Table closure used by the K-value path.
    ComponentNames componentNames{}; ///< Stable names used by diagnostics/output.
    int fullyCompositionalWaterComponent{-1}; ///< H2O component index for full three-phase flash seeding.
    PhaseNames phaseNames{}; ///< Phase labels in `Indices::Phase` ordering.

    // 外部算例接口：这里注入岩石–流体本构；饱和度为体积分数，
    // Natural 内核统一按 `lambda_alpha = k_r,alpha / mu_alpha` 组装相流度。
    // Explicit case-owned constitutive override. Input is phase MASS fraction.
    // Empty by default; existing LBC/IAPWS/McBride-Wright behavior is unchanged.
    std::function<ValueType(ValueType, const Composition &)> flowViscosityOverride;
    UnaryProperty gasRelativePermeability;   ///< `k_rg(S_g)`.
    UnaryProperty oilRelativePermeability;   ///< Two-phase/fallback `k_ro(S_o)`.
    UnaryProperty waterRelativePermeability; ///< `k_rw(S_w)`.
    UnaryProperty waterViscosity;            ///< `mu_w(p)` [Pa s].
    UnaryProperty waterFormationVolumeFactor; ///< `B_w(p)`, dimensionless.
    UnaryProperty poreVolumeMultiplier;      ///< Rock pore-volume multiplier `M_pv(p)`.
    ThreePhaseRelativePermeability threePhaseOilRelativePermeability; ///< `k_ro(S_o,S_g,S_w)`.

    int aqueousCO2Component{1}; ///< Component coupled to the optional aqueous-CO2 model.
    double temperature{298.15}; ///< Reservoir temperature, K.

private:
    template <class T, std::size_t N>
    static void copyList(
        std::initializer_list<T> source,
        std::array<T, N> &destination,
        const char *label)
    {
        if (source.size() != N)
            throw std::invalid_argument(std::string(label) + " list has wrong size.");
        std::copy(source.begin(), source.end(), destination.begin());
    }

    void validate() const
    {
        if (!(temperature > 0.0) || !std::isfinite(temperature))
            throw std::invalid_argument("FluidSystem temperature must be finite and positive.");
        for (double density : surfaceDensity)
            if (!(density > 0.0) || !std::isfinite(density))
                throw std::invalid_argument("Surface densities must be finite and positive.");
    }

    std::optional<AqueousCO2Model> aqueousCO2_;
    std::optional<CompetitiveLangmuirAdsorption<numComponents>> adsorption_;
    std::optional<Iapws2008IndustrialAqueousViscosity<numComponents>>
        iapws2008AqueousViscosity_;
    std::optional<McBrideWright2015AqueousViscosity<numComponents>>
        mcBrideWrightAqueousViscosity_;
};

} // namespace MPMC
