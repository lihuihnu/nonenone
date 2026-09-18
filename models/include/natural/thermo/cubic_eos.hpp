/**
 * @file cubic_eos.hpp
 * @brief PR/SRK/SW 立方状态方程的混合规则、根选择和逸度计算。
 */
#pragma once

#include <ad/Math.hpp>
#include <common/math.hpp>
#include <common/units.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/numerics.hpp>
#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/aqueous_volume.hpp>
#include <natural/thermo/cubic_plus_association.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <optional>

namespace MPMC
{

/**
 * @brief 通用两常数立方状态方程，默认用于 Peng-Robinson。
 *
 * 采用两常数立方 EOS
 * `p = RT/(v-b) - a alpha / [(v+m1 b)(v+m2 b)]`。
 * 无量纲混合参数使用
 * `A = sum_i sum_j x_i x_j A_ij`、`B = sum_i x_i B_i`，其中
 * `A_ij = sqrt(A_i A_j)(1-k_ij)`。组分逸度为 `f_i=phi_i x_i p`。
 * EOS 层只依赖组分数据与标量运算，不依赖 PETSc、网格或 AD 存储细节。
 *
 * 对 Peng-Robinson：
 *   omegaA = 0.4572355
 *   omegaB = 0.0779691
 *   m1 = 1 + sqrt(2)
 *   m2 = 1 - sqrt(2)
 */
template <class Indices>
class CubicEquationOfState
{
public:
    static constexpr int numComponents = Indices::numComponents;

    template <class Scalar>
    using Composition = std::array<Scalar, numComponents>;

    template <class Scalar>
    using InteractionMatrix =
        std::array<std::array<Scalar, numComponents>, numComponents>;

    template <class Scalar>
    struct MixingParameters
    {
        InteractionMatrix<Scalar> Aij{};
        Composition<Scalar> Bi{};
    };

    template <class Scalar>
    struct PhaseMixing
    {
        Composition<Scalar> Si{};
        Scalar A{0.0};
        Scalar B{0.0};
    };

    template <class Scalar>
    struct PhaseResult
    {
        Scalar compressibility{1.0};
        Composition<Scalar> fugacity{};
        Composition<Scalar> fugacityCoefficient{};
    };

    using ValueType = typename Indices::ValueType;
    using EquilibriumConstantFunction =
        std::function<ValueType(
            ValueType,
            ValueType,
            const Composition<ValueType> &)>;
    using EquilibriumConstantFunctions =
        std::array<EquilibriumConstantFunction, numComponents>;
    using BinaryInteractionFunction =
        std::function<double(int, int, double)>; // i, j, T [K]
    using AqueousWaterBipFunction =
        std::function<double(double, double)>; // T [K], NaCl molality [mol/kg H2O]
    using AqueousWaterBipFunctions =
        std::array<AqueousWaterBipFunction, numComponents>;

    struct SoreideWhitsonOptions
    {
        int waterComponent{-1};
        double salinityMolality{0.0};
        AqueousWaterBipFunctions aqueousWaterBip{};
        // Optional phase-specific volume correction [m3/mol].  It is applied
        // only to the thermodynamic Water role and therefore changes density
        // and phase volume without entering fugacity or phase-equilibrium
        // conditions.  A common use is a dissolved-CO2 correction fitted to
        // the Garcia apparent-molar-volume correlation.
        Composition<double> aqueousVolumeTranslation{};
    };

    /**
     * @brief SRK-CPA 热力学后端的纯组分和交叉缔合参数。
     *
     * 立方项使用 `a0`、`b`、`c1`；缔合位点用供体/受体数量表示。
     * 纯组分 `epsilon/beta` 默认按 CR-1 组合；显式交叉矩阵的正值会覆盖默认组合规则。
     */
    struct CubicPlusAssociationOptions
    {
        Composition<double> a0{};                 // Pa m6/mol2
        Composition<double> b{};                  // m3/mol
        Composition<double> c1{};                 // dimensionless
        Composition<double> associationEnergy{};  // J/mol
        Composition<double> associationVolume{};  // dimensionless
        std::array<int, numComponents> donorSites{};
        std::array<int, numComponents> acceptorSites{};
        InteractionMatrix<double> crossAssociationEnergy{};
        InteractionMatrix<double> crossAssociationVolume{};
        CpaCubicPhysicalTerm physicalTerm{
            CpaCubicPhysicalTerm::SoaveRedlichKwong};
        CpaRadialDistribution radialDistribution{CpaRadialDistribution::Simplified};
        // This bounds both the safeguarded site solver and its damped
        // fixed-point fallback. Converged states still exit early.
        int maximumAssociationIterations{2000};
        double associationTolerance{1.0e-12};
        double associationDamping{0.2};
    };

    struct AqueousVolumeOptions
    {
        int waterComponent{-1};
        int co2Component{-1};
        double waterMolarMass{0.0};
        double maximumUnsupportedMoleFraction{1.0e-4};
    };

    struct ThermodynamicProfile
    {
        std::uint64_t cpaMixingCalls{0};
        std::uint64_t cpaTemperatureCacheHits{0};
        std::uint64_t cpaAssociationCalls{0};
        std::uint64_t cpaAssociationAnalyticCalls{0};
        std::uint64_t cpaAssociationIterativeCalls{0};
        std::uint64_t cpaAssociationIterations{0};
        std::uint64_t cpaDensityRootCalls{0};
        std::uint64_t cpaDensityResidualEvaluations{0};
        std::uint64_t cpaDensityBisectionIterations{0};
        std::uint64_t cpaPhaseResultCalls{0};
        double cpaMixingSeconds{0.0};
        double cpaAssociationSeconds{0.0};
        double cpaDensityRootSeconds{0.0};
        double cpaPhaseResultSeconds{0.0};
        bool cpaTemperatureCacheEnabled{false};
        bool cpaWaterOnlyAnalyticFastPathAvailable{false};
    };

    CubicEquationOfState() = default;

    CubicEquationOfState(
        double omegaA,
        double omegaB,
        CompositionalMixture<Indices> mixture,
        int eosType,
        double m1,
        double m2,
        double minimumComposition = NaturalNumerics::minimumComposition)
        : omegaA_(omegaA),
          omegaB_(omegaB),
          mixture_(std::move(mixture)),
          eosType_(eosType),
          m1_(m1),
          m2_(m2),
          minimumComposition_(minimumComposition)
    {
        if (!(omegaA_ > 0.0) || !(omegaB_ > 0.0))
            throw std::invalid_argument("EOS omega constants must be positive.");
        if (m1_ == m2_)
            throw std::invalid_argument("Cubic EOS m1 and m2 must differ.");
        if (eosType_ != 1 && eosType_ != 5)
            throw std::invalid_argument("Natural CubicEquationOfState supports EOS types 1 and 5.");
    }

    /** @brief 返回 EOS 和黏度关联式共用的只读组分物性。 */
    [[nodiscard]] const CompositionalMixture<Indices> &mixture() const noexcept
    {
        return mixture_;
    }

    /** @brief 启用 Søreide–Whitson 的 PR 水/盐水修正。 */
    void configureSoreideWhitson(SoreideWhitsonOptions options)
    {
        if (options.waterComponent < 0 || options.waterComponent >= numComponents)
            throw std::invalid_argument("SW water-component index is out of range.");
        if (!(options.salinityMolality >= 0.0) ||
            !std::isfinite(options.salinityMolality))
        {
            throw std::invalid_argument("SW salinity must be finite and non-negative.");
        }
        for (double value : options.aqueousVolumeTranslation)
        {
            if (!std::isfinite(value))
                throw std::invalid_argument(
                    "SW aqueous volume-translation constants must be finite.");
        }
        aqueousVolumeTranslationEnabled_ = std::any_of(
            options.aqueousVolumeTranslation.begin(),
            options.aqueousVolumeTranslation.end(),
            [](double value) { return value != 0.0; });
        thermodynamicModel_ = CubicThermodynamicModel::SoreideWhitson;
        soreideWhitson_ = std::move(options);
        cpaTemperatureCache_.reset();
    }

    /** @brief 启用 SRK-CPA 热力学后端。 */
    void configureCubicPlusAssociation(CubicPlusAssociationOptions options)
    {
        if (options.maximumAssociationIterations <= 0)
            throw std::invalid_argument("CPA association iteration limit must be positive.");
        if (!(options.associationTolerance > 0.0) ||
            !std::isfinite(options.associationTolerance))
            throw std::invalid_argument("CPA association tolerance must be finite and positive.");
        if (!(options.associationDamping > 0.0) ||
            !(options.associationDamping <= 1.0) ||
            !std::isfinite(options.associationDamping))
            throw std::invalid_argument("CPA association damping must be in (0,1].");

        // Validate all explicit cross-association entries before deciding
        // whether a site-bearing component is self-associating or
        // cross-association-only (solvating).  CPA formulations for aromatic
        // hydrocarbons commonly use the latter: the inert component has
        // epsilon_AiBi = 0 but carries a site that associates explicitly with
        // a donor/acceptor on water.
        for (int i = 0; i < numComponents; ++i)
        {
            for (int j = 0; j < numComponents; ++j)
            {
                const auto ii = static_cast<std::size_t>(i);
                const auto jj = static_cast<std::size_t>(j);
                const double epsilon =
                    options.crossAssociationEnergy[ii][jj];
                const double beta =
                    options.crossAssociationVolume[ii][jj];
                if (!std::isfinite(epsilon) || !std::isfinite(beta) ||
                    epsilon < 0.0 || beta < 0.0)
                    throw std::invalid_argument(
                        "CPA explicit cross-association parameters must be finite and non-negative.");
                if ((epsilon > 0.0) != (beta > 0.0))
                    throw std::invalid_argument(
                        "CPA explicit cross-association requires epsilon and beta together.");
            }
        }

        for (int i = 0; i < numComponents; ++i)
        {
            const auto idx = static_cast<std::size_t>(i);
            if (!(options.a0[idx] > 0.0) || !(options.b[idx] > 0.0) ||
                !std::isfinite(options.a0[idx]) || !std::isfinite(options.b[idx]) ||
                !std::isfinite(options.c1[idx]))
                throw std::invalid_argument("CPA cubic pure-component parameters must be finite and positive where required.");
            if (options.donorSites[idx] < 0 || options.acceptorSites[idx] < 0)
                throw std::invalid_argument("CPA association site counts cannot be negative.");
            if (!std::isfinite(options.associationEnergy[idx]) ||
                !std::isfinite(options.associationVolume[idx]) ||
                options.associationEnergy[idx] < 0.0 ||
                options.associationVolume[idx] < 0.0)
                throw std::invalid_argument(
                    "CPA pure association epsilon and beta must be finite and non-negative.");

            const bool hasSites =
                options.donorSites[idx] + options.acceptorSites[idx] > 0;
            if (!hasSites)
                continue;

            const bool hasPureEnergy = options.associationEnergy[idx] > 0.0;
            const bool hasPureVolume = options.associationVolume[idx] > 0.0;
            if (hasPureEnergy != hasPureVolume)
                throw std::invalid_argument(
                    "CPA self-association requires pure epsilon and beta together.");
            if (hasPureEnergy)
                continue;

            // Cross-association-only site carriers are valid only when every
            // declared site direction has at least one explicit positive
            // partner.  This prevents a zero pure epsilon from silently
            // disabling a declared site while allowing faithful solvating
            // mixtures such as water + aromatic pseudo-components.
            bool donorCovered = options.donorSites[idx] == 0;
            bool acceptorCovered = options.acceptorSites[idx] == 0;
            for (int j = 0; j < numComponents; ++j)
            {
                const auto jj = static_cast<std::size_t>(j);
                if (!donorCovered && options.acceptorSites[jj] > 0 &&
                    options.crossAssociationEnergy[idx][jj] > 0.0 &&
                    options.crossAssociationVolume[idx][jj] > 0.0)
                    donorCovered = true;
                if (!acceptorCovered && options.donorSites[jj] > 0 &&
                    options.crossAssociationEnergy[jj][idx] > 0.0 &&
                    options.crossAssociationVolume[jj][idx] > 0.0)
                    acceptorCovered = true;
            }
            if (!donorCovered || !acceptorCovered)
                throw std::invalid_argument(
                    "CPA cross-association-only component has an uncovered declared site.");
        }

        thermodynamicModel_ = CubicThermodynamicModel::CubicPlusAssociation;
        soreideWhitson_.reset();
        cubicPlusAssociation_ = std::move(options);
        cpaTemperatureCache_.reset();
        resetThermodynamicProfile();
    }

    /** @brief 恢复普通 Peng–Robinson 热力学后端。 */
    void configurePengRobinson() noexcept
    {
        thermodynamicModel_ = CubicThermodynamicModel::PengRobinson;
        soreideWhitson_.reset();
        cubicPlusAssociation_.reset();
        cpaTemperatureCache_.reset();
    }

    [[nodiscard]] CubicThermodynamicModel thermodynamicModel() const noexcept
    {
        return thermodynamicModel_;
    }

    [[nodiscard]] bool usesSoreideWhitson() const noexcept
    {
        return thermodynamicModel_ == CubicThermodynamicModel::SoreideWhitson;
    }

    [[nodiscard]] bool usesCubicPlusAssociation() const noexcept
    {
        return thermodynamicModel_ == CubicThermodynamicModel::CubicPlusAssociation;
    }


    /** @brief 启用或关闭低开销热力学性能统计。 */
    void setThermodynamicProfilerEnabled(bool enabled) const noexcept
    {
        thermodynamicProfilerEnabled_ = enabled;
    }

    [[nodiscard]] bool thermodynamicProfilerEnabled() const noexcept
    {
        return thermodynamicProfilerEnabled_;
    }

    void resetThermodynamicProfile() const noexcept
    {
        thermodynamicProfile_ = ThermodynamicProfile{};
        thermodynamicProfile_.cpaTemperatureCacheEnabled = cpaTemperatureCache_.has_value();
        thermodynamicProfile_.cpaWaterOnlyAnalyticFastPathAvailable =
            cpaWaterOnlyAnalyticFastPathAvailable_();
    }

    [[nodiscard]] ThermodynamicProfile thermodynamicProfile() const noexcept
    {
        auto profile = thermodynamicProfile_;
        profile.cpaTemperatureCacheEnabled = cpaTemperatureCache_.has_value();
        profile.cpaWaterOnlyAnalyticFastPathAvailable =
            cpaWaterOnlyAnalyticFastPathAvailable_();
        return profile;
    }

    /**
     * @brief 为等温 CPA 运行缓存只依赖温度的纯组分和交叉参数。
     *
     * 储层模型没有温度主未知量，因此这里只缓存依赖 T 的纯组分/交叉参数，
     * 所有组成相关 mixing 仍保持可微并在调用时实时计算。
     */
    void configureCpaTemperatureCache(double temperature)
    {
        if (!usesCubicPlusAssociation())
            throw std::logic_error("CPA temperature cache requires CPA configuration.");
        if (!(temperature > 0.0) || !std::isfinite(temperature))
            throw std::invalid_argument("CPA cache temperature must be finite and positive.");
        cpaTemperatureCache_ = buildCpaTemperatureCache_(temperature);
        thermodynamicProfile_.cpaTemperatureCacheEnabled = true;
    }

    void clearCpaTemperatureCache() noexcept
    {
        cpaTemperatureCache_.reset();
        thermodynamicProfile_.cpaTemperatureCacheEnabled = false;
    }

    /**
     * @brief 配置可选的温度相关立方 EOS 二元作用系数 `k_ij(T)`。
     *
     * 该回调用于普通 PR、SW 的非水相部分以及 CPA 的 SRK 立方项；未配置时使用常数 BIP 矩阵。
     * 当前储层求解等温，但离线 PVT/相图工具允许扫描温度，因此接口保留 `T`。
     */
    void configureBinaryInteractionFunction(BinaryInteractionFunction function)
    {
        binaryInteractionFunction_ = std::move(function);
        cpaTemperatureCache_.reset();
    }

    void clearBinaryInteractionFunction() noexcept
    {
        binaryInteractionFunction_ = {};
        cpaTemperatureCache_.reset();
    }

    [[nodiscard]] bool usesTemperatureDependentBinaryInteraction() const noexcept
    {
        return static_cast<bool>(binaryInteractionFunction_);
    }

    /** @brief 返回温度 `T` [K] 下实际使用的立方 EOS `k_ij`。 */
    [[nodiscard]] double binaryInteractionCoefficient(
        int i, int j, double temperature) const
    {
        const double kij = binaryInteractionFunction_
            ? binaryInteractionFunction_(i, j, temperature)
            : mixture_.binaryInteractionCoefficient(i, j);
        if (!std::isfinite(kij))
            throw std::runtime_error("Cubic EOS BIP correlation returned a non-finite value.");

        // 当前 classical one-fluid mixing rule 及其逸度导数使用 A_ij=A_ji。
        // 自定义温度相关回调若破坏 k_ij=k_ji，会使 A 的组成导数与实现中的 2*S_i
        // 不一致，因此在首次半矩阵访问时明确拒绝，而不是静默产生错误化学势。
        if (binaryInteractionFunction_ && i < j)
        {
            const double kji = binaryInteractionFunction_(j, i, temperature);
            if (!std::isfinite(kji))
                throw std::runtime_error("Cubic EOS BIP correlation returned a non-finite reverse value.");
            const double scale = std::max({1.0, std::abs(kij), std::abs(kji)});
            if (std::abs(kij - kji) > 1.0e-12 * scale)
                throw std::invalid_argument(
                    "Cubic EOS temperature-dependent BIP must satisfy k_ij(T)=k_ji(T)." );
        }
        return kij;
    }

    /**
     * @brief 配置组分体积平移常数 `c_i` [m3/mol]。
     *
     * 数学：`v = ZRT/P - sum_i x_i c_i`。体积平移只修正密度/相体积，
     * 不改变 Z、逸度系数和相平衡条件；全零数组严格对应旧行为。
     */
    void configureVolumeTranslation(Composition<double> componentShift)
    {
        for (double value : componentShift)
        {
            if (!std::isfinite(value))
                throw std::invalid_argument("Volume-translation constants must be finite.");
        }
        volumeTranslationEnabled_ = std::any_of(
            componentShift.begin(), componentShift.end(),
            [](double value) { return value != 0.0; });
        volumeTranslation_ = std::move(componentShift);
    }

    /** @brief 仅对 Water 相角色启用 IAPWS+Garcia 摩尔体积闭包。 */
    void configureAqueousVolume(AqueousVolumeOptions options)
    {
        aqueousVolumeModel_.emplace(
            options.waterComponent,
            options.co2Component,
            options.waterMolarMass,
            options.maximumUnsupportedMoleFraction);
    }

    /**
     * @brief 将公共 Water 角色限制在已验证的水相物性组成域内。
     *
     * A hydrocarbon-rich phase can contain more H2O on a mole basis than any
     * individual heavy component.  Mole-fraction dominance alone must not move
     * that phase into the Water slot, because doing so also selects aqueous
     * density/viscosity closures.  This independent domain guard is useful when
     * an aqueous transport property is configured without a special volume
     * model.
     */
    void configureAqueousCompositionDomain(
        int waterComponent,
        double maximumSoluteMoleFraction)
    {
        if (waterComponent < 0 || waterComponent >= numComponents ||
            !(maximumSoluteMoleFraction >= 0.0) ||
            !std::isfinite(maximumSoluteMoleFraction))
        {
            throw std::invalid_argument(
                "Invalid aqueous composition-domain configuration.");
        }
        aqueousDomainWaterComponent_ = waterComponent;
        aqueousDomainMaximumSoluteMoleFraction_ = maximumSoluteMoleFraction;
    }

    [[nodiscard]] bool usesAqueousVolume() const noexcept
    {
        return aqueousVolumeModel_.has_value();
    }

    /**
     * @brief 检查 Water-role 组成是否满足全部已配置水相物性适用域。
     *
     * 未配置专用水相体积闭包时所有有限 EOS 组成均由原 cubic 路径处理。
     */
    [[nodiscard]] bool aqueousVolumeCompositionSupported(
        const Composition<double> &composition) const noexcept
    {
        if (aqueousDomainWaterComponent_ >= 0)
        {
            double solute = 0.0;
            for (int component = 0; component < numComponents; ++component)
            {
                if (component != aqueousDomainWaterComponent_)
                    solute += composition[static_cast<std::size_t>(component)];
            }
            if (!std::isfinite(solute) ||
                solute > aqueousDomainMaximumSoluteMoleFraction_)
            {
                return false;
            }
        }
        return !aqueousVolumeModel_ ||
            aqueousVolumeModel_->compositionSupported(composition);
    }

    [[nodiscard]] const Composition<double> &volumeTranslation() const noexcept
    {
        return volumeTranslation_;
    }

    [[nodiscard]] bool usesVolumeTranslation() const noexcept
    {
        return volumeTranslationEnabled_ || aqueousVolumeTranslationEnabled_ ||
            usesAqueousVolume();
    }

    template <class Scalar>
    [[nodiscard]] Scalar translatedMolarVolume(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        const Scalar &compressibility,
        CompositionalPhase phase) const
    {
        constexpr double R = units::gasConstant;
        if (phase == CompositionalPhase::Water && aqueousVolumeModel_)
            return aqueousVolumeModel_->molarVolume(
                pressure, temperature, composition);
        Scalar shift = 0.0;
        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t component = static_cast<std::size_t>(i);
            double componentShift = volumeTranslation_[component];
            if (phase == CompositionalPhase::Water && usesSoreideWhitson())
                componentShift += soreideWhitsonOptions_()
                    .aqueousVolumeTranslation[component];
            shift += composition[component] * componentShift;
        }
        const Scalar volume = compressibility * R * temperature / pressure - shift;
        if (!(scalarValue(volume) > 0.0) || !std::isfinite(scalarValue(volume)))
            throw std::runtime_error("Volume translation produced a non-positive phase molar volume.");
        return volume;
    }

    /** @brief 兼容旧接口：无相角色时只应用通用组分体积平移。 */
    template <class Scalar>
    [[nodiscard]] Scalar translatedMolarVolume(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        const Scalar &compressibility) const
    {
        return translatedMolarVolume(
            pressure, temperature, composition, compressibility,
            CompositionalPhase::Oil);
    }

    template <class Scalar>
    [[nodiscard]] Scalar molarDensity(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        const Scalar &compressibility,
        CompositionalPhase phase) const
    {
        constexpr double R = units::gasConstant;
        if (!usesVolumeTranslation())
            return pressure / (R * temperature * compressibility);
        return Scalar(1.0) / translatedMolarVolume(
            pressure, temperature, composition, compressibility, phase);
    }

    template <class Scalar>
    [[nodiscard]] Scalar molarDensity(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        const Scalar &compressibility) const
    {
        constexpr double R = units::gasConstant;
        if (!usesVolumeTranslation())
            return pressure / (R * temperature * compressibility);
        return Scalar(1.0) / translatedMolarVolume(
            pressure, temperature, composition, compressibility,
            CompositionalPhase::Oil);
    }

    /** @brief 返回只读 CPA 参数，供离线诊断和回归工具使用。 */
    [[nodiscard]] const CubicPlusAssociationOptions &cubicPlusAssociationOptions() const
    {
        return cpaOptions_();
    }

    /**
     * @brief 为旧 K-value 模式配置每个组分的 `K_i(p,T,z)` 关联式。
     *
     * 配置后相平衡不再使用立方 EOS 逸度相等；定长数组保证全部组分同时配置，避免半配置状态。
     */
    void configureEquilibriumConstants(EquilibriumConstantFunctions functions)
    {
        equilibriumConstantFunctions_ = std::move(functions);
    }

    /** @brief 判断当前是否使用 K-value 模式而非 EOS 逸度 flash。 */
    [[nodiscard]] bool usesEquilibriumConstants() const noexcept
    {
        return equilibriumConstantFunctions_.has_value();
    }

    /**
     * @brief 构造指定热力学相在 `(p,T)` 下的无量纲 `A_i/B_i/A_ij`。
     *
     * 物理：普通 PR 各相共享同一 BIP 和 PR alpha；SW 仍使用 PR 立方常数，
     * 但所有相采用 SW 的 H2O alpha，并且**只有 Aqueous 角色**使用水相 H2O–溶质 BIP。
     * 因此 SW 调用方必须在非线性迭代期间保持热力学角色固定。
     */
    template <class Scalar>
    [[nodiscard]] MixingParameters<Scalar>
    mixingParameters(
        const Scalar &pressure,
        const Scalar &temperature,
        CompositionalPhase phase) const
    {
        MixingParameters<Scalar> result;
        Composition<Scalar> sqrtAi{};

        for (int i = 0; i < numComponents; ++i)
        {
            const Scalar pr = pressure / mixture_.criticalPressure(i);
            const Scalar tr = temperature / mixture_.criticalTemperature(i);
            Scalar alpha = 1.0;

            if (usesSoreideWhitson() &&
                i == soreideWhitsonOptions_().waterComponent)
            {
                alpha = SoreideWhitsonCorrelations::waterAlpha(
                    temperature, mixture_.criticalTemperature(i),
                    soreideWhitsonOptions_().salinityMolality);
            }
            else
            {
                const double omega = mixture_.acentricFactor(i);
                double m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
                if (eosType_ == 5 && omega > 0.49)
                {
                    m = 0.379642 + 1.48503 * omega
                        - 0.164423 * omega * omega
                        + 0.016666 * omega * omega * omega;
                }
                const Scalar alphaRoot =
                    Scalar(1.0) + m * (Scalar(1.0) - sqrt(tr));
                alpha = alphaRoot * alphaRoot;
            }

            const Scalar reducedA = omegaA_ * alpha * pr / (tr * tr);
            sqrtAi[static_cast<std::size_t>(i)] = sqrt(reducedA);
            result.Bi[static_cast<std::size_t>(i)] = omegaB_ * pr / tr;
        }

        const double temperatureValue = scalarValue(temperature);
        const bool useAqueousSwBip =
            usesSoreideWhitson() && phase == CompositionalPhase::Water;

        // Classical one-fluid mixing requires k_ij = k_ji, so A_ij is
        // symmetric. Build the upper triangle once and mirror it to avoid
        // duplicate A_ij construction. Temperature-dependent callbacks still
        // evaluate the reverse entry when needed to enforce the symmetry contract.
        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t ii = static_cast<std::size_t>(i);
            for (int j = i; j < numComponents; ++j)
            {
                const std::size_t jj = static_cast<std::size_t>(j);
                double kij = binaryInteractionCoefficient(i, j, temperatureValue);
                // 物理：SW 仅在固定 Aqueous 热力学角色上替换 H2O–溶质 BIP；
                // 这里不根据瞬时含水量判断，避免 EOS 参数随非线性迭代发生离散跳变。
                if (useAqueousSwBip)
                    kij = soreideWhitsonAqueousBip_(i, j, temperatureValue);

                const Scalar aij = sqrtAi[ii] * sqrtAi[jj] * (1.0 - kij);
                result.Aij[ii][jj] = aij;
                result.Aij[jj][ii] = aij;
            }
        }
        return result;
    }

    /** @brief 兼容旧接口：按普通非水液相角色构造混合参数。 */
    template <class Scalar>
    [[nodiscard]] MixingParameters<Scalar>
    mixingParameters(const Scalar &pressure, const Scalar &temperature) const
    {
        return mixingParameters(pressure, temperature, CompositionalPhase::Oil);
    }

    /**
     * @brief 按给定相组成混合组分 EOS 参数。
     *
     * 数学：`A=sum_i sum_j x_i x_j A_ij`，`B=sum_i x_i B_i`，
     * `S_i=sum_j x_j A_ij`。
     */
    template <class Scalar>
    [[nodiscard]] PhaseMixing<Scalar>
    phaseMixing(
        const Composition<Scalar> &composition,
        const MixingParameters<Scalar> &parameters) const
    {
        PhaseMixing<Scalar> result;
        for (int i = 0; i < numComponents; ++i)
        {
            Scalar si = 0.0;
            for (int j = 0; j < numComponents; ++j)
            {
                si += composition[static_cast<std::size_t>(j)] *
                      parameters.Aij[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            }
            result.Si[static_cast<std::size_t>(i)] = si;
            result.B += composition[static_cast<std::size_t>(i)] *
                        parameters.Bi[static_cast<std::size_t>(i)];
        }

        for (int i = 0; i < numComponents; ++i)
        {
            result.A += composition[static_cast<std::size_t>(i)] *
                        result.Si[static_cast<std::size_t>(i)];
        }
        return result;
    }

private:
    template <class Scalar>
    [[nodiscard]] std::array<Scalar, 3>
    cubicCoefficients(const Scalar &A, const Scalar &B) const
    {
        const Scalar e0 = -(A * B + m1_ * m2_ * B * B * (B + 1.0));
        const Scalar e1 = A
            - (m1_ + m2_ - m1_ * m2_) * B * B
            - (m1_ + m2_) * B;
        const Scalar e2 = (m1_ + m2_ - 1.0) * B - 1.0;
        return {e2, e1, e0};
    }

    template <class Scalar>
    [[nodiscard]] std::array<Scalar, 3>
    solveCubic(const Scalar &A, const Scalar &B) const
    {
        const auto coefficients = cubicCoefficients(A, B);
        const Scalar a = coefficients[0];
        const Scalar b = coefficients[1];
        const Scalar c = coefficients[2];
        const Scalar Q = (a * a - 3.0 * b) / 9.0;
        const Scalar R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;
        const Scalar discriminant = R * R - Q * Q * Q;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        std::array<Scalar, 3> roots{Scalar(nan), Scalar(nan), Scalar(nan)};

        if (scalarValue(discriminant) < 0.0)
        {
            constexpr double pi = 3.141592653589793238462643383279502884;
            Scalar acosArgument = R / sqrt(Q * Q * Q);
            // 临界/重根附近理论值应位于 [-1,1]，浮点舍入可能产生极小越界。
            // 只在越界时投影，避免 acos 返回 NaN；正常区间内保留完整 AD 导数。
            if (scalarValue(acosArgument) > 1.0)
                acosArgument = Scalar(1.0);
            else if (scalarValue(acosArgument) < -1.0)
                acosArgument = Scalar(-1.0);
            const Scalar theta = acos(acosArgument);
            roots[0] = -2.0 * sqrt(Q) * cos(theta / 3.0) - a / 3.0;
            roots[1] = -2.0 * sqrt(Q) * cos((theta + 2.0 * pi) / 3.0) - a / 3.0;
            roots[2] = -2.0 * sqrt(Q) * cos((theta - 2.0 * pi) / 3.0) - a / 3.0;
        }
        else
        {
            Scalar S;
            if (scalarValue(R) > 0.0)
                S = -pow(abs(R) + sqrt(discriminant), 1.0 / 3.0);
            else
                S = pow(abs(R) + sqrt(discriminant), 1.0 / 3.0);

            Scalar T = 0.0;
            if (std::abs(scalarValue(S)) > 0.0)
                T = Q / S;
            roots[0] = S + T - a / 3.0;
        }

        for (auto &root : roots)
        {
            if (std::isfinite(scalarValue(root)) && scalarValue(root) <= 0.0)
                root = Scalar(nan);
        }
        return roots;
    }

public:
    /** @brief 返回满足 `Z>B` 的最小可接受根，作为液相型立方 EOS 根。 */
    template <class Scalar>
    [[nodiscard]] Scalar liquidRoot(const Scalar &A, const Scalar &B) const
    {
        const auto roots = solveCubic(A, B);
        bool found = false;
        Scalar best = std::numeric_limits<double>::quiet_NaN();
        for (const auto &root : roots)
        {
            const double z = scalarValue(root);
            if (!std::isfinite(z) || !(z > scalarValue(B)))
                continue;
            if (!found || z < scalarValue(best))
            {
                best = root;
                found = true;
            }
        }
        if (!found)
            throw std::runtime_error("Cubic EOS has no admissible liquid compressibility root.");
        return best;
    }

    /** @brief 返回满足 `Z>B` 的最大可接受根，作为气相型立方 EOS 根。 */
    template <class Scalar>
    [[nodiscard]] Scalar vaporRoot(const Scalar &A, const Scalar &B) const
    {
        const auto roots = solveCubic(A, B);
        bool found = false;
        Scalar best = std::numeric_limits<double>::quiet_NaN();
        for (const auto &root : roots)
        {
            const double z = scalarValue(root);
            if (!std::isfinite(z) || !(z > scalarValue(B)))
                continue;
            if (!found || z > scalarValue(best))
            {
                best = root;
                found = true;
            }
        }
        if (!found)
            throw std::runtime_error("Cubic EOS has no admissible vapor compressibility root.");
        return best;
    }

    /** @brief 计算单相的压缩因子、逸度系数和逸度。 */
    template <class Scalar>
    [[nodiscard]] PhaseResult<Scalar>
    phaseResult(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        bool liquid,
        bool selectGibbsMinimum = false) const
    {
        if (usesCubicPlusAssociation())
            return cpaPhaseResult_(pressure, temperature, composition, liquid, selectGibbsMinimum);

        const auto parameters = mixingParameters(
            pressure,
            temperature,
            liquid ? CompositionalPhase::Oil : CompositionalPhase::Gas);
        const auto mix = phaseMixing(composition, parameters);
        Scalar Z = liquid ? liquidRoot(mix.A, mix.B) : vaporRoot(mix.A, mix.B);

        if (selectGibbsMinimum)
            Z = gibbsMinimumRoot(pressure, composition, parameters, mix);

        auto result = fugacityAtRoot(pressure, composition, parameters, mix, Z);
        result.compressibility = Z;
        return result;
    }

    /** @brief 全组分 O/G/W 路径使用的“相角色感知”EOS 评价接口。 */
    template <class Scalar>
    [[nodiscard]] PhaseResult<Scalar>
    phaseResult(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        CompositionalPhase phase,
        bool selectGibbsMinimum = false) const
    {
        if (usesCubicPlusAssociation())
            return cpaPhaseResult_(
                pressure, temperature, composition,
                phase != CompositionalPhase::Gas, selectGibbsMinimum);

        const auto parameters = mixingParameters(pressure, temperature, phase);
        const auto mix = phaseMixing(composition, parameters);
        const bool liquid = phase != CompositionalPhase::Gas;
        Scalar Z = liquid ? liquidRoot(mix.A, mix.B) : vaporRoot(mix.A, mix.B);
        if (selectGibbsMinimum)
            Z = gibbsMinimumRoot(pressure, composition, parameters, mix);
        auto result = fugacityAtRoot(pressure, composition, parameters, mix, Z);
        result.compressibility = Z;
        return result;
    }

    /** @brief 只计算选定的压缩因子根，不计算逐组分逸度。
     *
     * 在不要求 Gibbs 最小根时与 `phaseResult(...).compressibility` 数学等价；
     * secondary-state 更新只需要 Z，因此该接口避免提前做随后会丢弃的逐组分逸度计算。
     */
    template <class Scalar>
    [[nodiscard]] Scalar compressibility(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        CompositionalPhase phase) const
    {
        const bool liquid = phase != CompositionalPhase::Gas;
        if (usesCubicPlusAssociation())
        {
            const Scalar rho = cpaDensityRoot_(
                pressure, temperature, composition, liquid);
            return pressure / (rho * cpaGasConstant_ * temperature);
        }

        const auto parameters = mixingParameters(pressure, temperature, phase);
        const auto mix = phaseMixing(composition, parameters);
        return liquid ? liquidRoot(mix.A, mix.B) : vaporRoot(mix.A, mix.B);
    }

    /** @brief 兼容旧油/气布尔参数的压缩因子接口。 */
    template <class Scalar>
    [[nodiscard]] Scalar compressibility(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        bool liquid) const
    {
        return compressibility(
            pressure,
            temperature,
            composition,
            liquid ? CompositionalPhase::Oil : CompositionalPhase::Gas);
    }

    /**
     * @brief 在给定压缩因子根 Z 上计算 `phi_i` 与 `f_i = phi_i x_i p`。
     */
    template <class Scalar>
    [[nodiscard]] PhaseResult<Scalar>
    fugacityAtRoot(
        const Scalar &pressure,
        const Composition<Scalar> &composition,
        const MixingParameters<Scalar> &parameters,
        const PhaseMixing<Scalar> &mix,
        const Scalar &Z) const
    {
        if (!(scalarValue(mix.B) > 0.0) || !(scalarValue(mix.A) > 0.0))
            throw std::runtime_error("Cubic EOS fugacity requires positive A and B.");
        if (!(scalarValue(Z - mix.B) > 0.0))
            throw std::runtime_error("Cubic EOS fugacity requires Z > B.");

        PhaseResult<Scalar> result;
        result.compressibility = Z;

        const Scalar commonA = -log(Z - mix.B);
        const Scalar logRatio = log((Z + m2_ * mix.B) / (Z + m1_ * mix.B));
        const Scalar commonB = logRatio * mix.A / ((m1_ - m2_) * mix.B);
        const Scalar commonC = (Z - 1.0) / mix.B;

        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            const Scalar lnPhi = commonA
                + commonB * (2.0 * mix.Si[idx] / mix.A - parameters.Bi[idx] / mix.B)
                + parameters.Bi[idx] * commonC;
            result.fugacityCoefficient[idx] = exp(lnPhi);
            result.fugacity[idx] = result.fugacityCoefficient[idx] * pressure * composition[idx];
        }
        return result;
    }

    /** @brief 在给定摩尔密度 [mol/m3] 下计算 CPA 压力，供诊断/回归使用。 */
    template <class Scalar>
    [[nodiscard]] Scalar cpaPressureAtMolarDensity(
        const Scalar &temperature,
        const Scalar &molarDensity,
        const Composition<Scalar> &composition) const
    {
        if (!usesCubicPlusAssociation())
            throw std::logic_error("CPA density diagnostic requires the CPA backend.");
        return cpaPressureAtDensity_(molarDensity, temperature, composition);
    }

    /** @brief 在固定密度下返回 CPA 供体/受体未缔合位点分数。 */
    template <class Scalar>
    [[nodiscard]] std::pair<Composition<Scalar>, Composition<Scalar>>
    cpaSiteFractionsAtMolarDensity(
        const Scalar &temperature,
        const Scalar &molarDensity,
        const Composition<Scalar> &composition) const
    {
        if (!usesCubicPlusAssociation())
            throw std::logic_error("CPA site diagnostic requires the CPA backend.");
        const auto mix = cpaMixing_(temperature, composition);
        const auto association = cpaAssociationState_(
            molarDensity, composition, mix);
        return {association.donorX, association.acceptorX};
    }

    /** @brief 按组分摩尔质量将相内摩尔分数换算为质量分数。 */
    template <class Scalar>
    [[nodiscard]] Composition<Scalar>
    massFractions(const Composition<Scalar> &moleFractions) const
    {
        Composition<Scalar> result{};
        Scalar totalMass = 0.0;
        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            result[idx] = moleFractions[idx] * mixture_.molecularWeight(i);
            totalMass += result[idx];
        }
        if (!(scalarValue(totalMass) > 0.0))
            throw std::runtime_error("Mass-fraction conversion requires positive mixture molar mass.");
        for (auto &entry : result)
            entry /= totalMass;
        return result;
    }

private:
    /**
     * @brief 求解 Rachford-Rice 中的 vapor fraction V=1-L。
     *
     * 使用有界 Newton，并在 Newton 步离开区间时回退到二分法。
     * 区间始终保持在物理范围 `[0,1]`。
     */
    template <class Scalar>
    [[nodiscard]] Scalar solveVaporFraction(
        Scalar initialLiquidFraction,
        const Composition<Scalar> &K,
        Composition<Scalar> z) const
    {
        normalizeComposition(z, minimumComposition_);

        double lower = 0.0;
        double upper = 1.0;
        Scalar V = 1.0 - initialLiquidFraction;
        double vValue = std::clamp(scalarValue(V), lower, upper);
        V = V + (vValue - scalarValue(V));

        auto residual = [&](const Scalar &v) {
            Scalar f = 0.0;
            Scalar positiveDenominator = 0.0;
            for (int i = 0; i < numComponents; ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(i);
                const Scalar km1 = K[idx] - 1.0;
                const Scalar denom = 1.0 + v * km1;
                f += z[idx] * km1 / denom;
                positiveDenominator += z[idx] * km1 * km1 / (denom * denom);
            }
            return std::pair<Scalar, Scalar>{f, positiveDenominator};
        };

        const Scalar fLower = residual(Scalar(lower)).first;
        const Scalar fUpper = residual(Scalar(upper)).first;
        if (scalarValue(fLower) <= 0.0)
            return Scalar(0.0);
        if (scalarValue(fUpper) >= 0.0)
            return Scalar(1.0);

        for (int iteration = 0; iteration < 100; ++iteration)
        {
            auto [f, denom] = residual(V);
            const double fv = scalarValue(f);
            if (std::abs(fv) <= 1.0e-12)
                return V;

            if (fv > 0.0)
                lower = scalarValue(V);
            else
                upper = scalarValue(V);

            Scalar candidate = V;
            if (scalarValue(denom) > 0.0)
                candidate = V + f / denom; // f'=-denom

            const double cv = scalarValue(candidate);
            if (!std::isfinite(cv) || cv <= lower || cv >= upper)
            {
                const double midpoint = 0.5 * (lower + upper);
                candidate = V + (midpoint - scalarValue(V));
            }
            V = candidate;
        }
        throw std::runtime_error("Rachford-Rice did not converge within 100 safeguarded iterations.");
    }

public:
    /**
     * @brief 求解两相 Rachford–Rice 方程并返回液相分率 L。
     *
     * 数学：`K_i=y_i/x_i`、`V=1-L`，求解
     * `sum_i z_i (K_i-1)/(1+V(K_i-1)) = 0`。
     */
    template <class Scalar>
    [[nodiscard]] Scalar solveRachfordRice(
        const Scalar &liquidFraction,
        const Composition<Scalar> &K,
        const Composition<Scalar> &z) const
    {
        return 1.0 - solveVaporFraction(liquidFraction, K, z);
    }

    /** @brief 由 `x_i=z_i/[L+(1-L)K_i]` 计算并归一化液相组成。 */
    template <class Scalar>
    [[nodiscard]] Composition<Scalar>
    liquidComposition(
        const Scalar &liquidFraction,
        const Composition<Scalar> &K,
        const Composition<Scalar> &z) const
    {
        Composition<Scalar> x{};
        Scalar sum = 0.0;
        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            const Scalar denom = liquidFraction + (1.0 - liquidFraction) * K[idx];
            x[idx] = z[idx] / denom;
            sum += x[idx];
        }
        if (!(scalarValue(sum) > 0.0))
            throw std::runtime_error("Liquid flash composition has zero total.");
        for (auto &entry : x)
            entry /= sum;
        return x;
    }

    /** @brief 由 `y_i=K_i z_i/[L+(1-L)K_i]` 计算并归一化气相组成。 */
    template <class Scalar>
    [[nodiscard]] Composition<Scalar>
    vaporComposition(
        const Scalar &liquidFraction,
        const Composition<Scalar> &K,
        const Composition<Scalar> &z) const
    {
        Composition<Scalar> y{};
        Scalar sum = 0.0;
        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            const Scalar denom = liquidFraction + (1.0 - liquidFraction) * K[idx];
            y[idx] = K[idx] * z[idx] / denom;
            sum += y[idx];
        }
        if (!(scalarValue(sum) > 0.0))
            throw std::runtime_error("Vapor flash composition has zero total.");
        for (auto &entry : y)
            entry /= sum;
        return y;
    }

    /** @brief 计算当前标量类型下全部已配置的 `K_i(p,T,z)`。 */
    template <class Scalar>
    [[nodiscard]] Composition<Scalar>
    evaluateEquilibriumConstants(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &z) const
    {
        if (!equilibriumConstantFunctions_)
            throw std::logic_error("Equilibrium-constant correlations are not configured.");

        Composition<Scalar> result{};
        Composition<ValueType> composition{};
        for (int i = 0; i < numComponents; ++i)
            composition[static_cast<std::size_t>(i)] =
                ValueType(z[static_cast<std::size_t>(i)]);

        for (int i = 0; i < numComponents; ++i)
        {
            const auto value = (*equilibriumConstantFunctions_)[static_cast<std::size_t>(i)](
                ValueType(pressure), ValueType(temperature), composition);
            result[static_cast<std::size_t>(i)] = MPMC::decay<Scalar>(value);
        }
        return result;
    }

private:
    template <class Scalar>
    [[nodiscard]] Scalar gibbsMinimumRoot(
        const Scalar &pressure,
        const Composition<Scalar> &composition,
        const MixingParameters<Scalar> &parameters,
        const PhaseMixing<Scalar> &mix) const
    {
        const auto roots = solveCubic(mix.A, mix.B);
        bool found = false;
        Scalar best = std::numeric_limits<double>::quiet_NaN();
        double bestGibbs = std::numeric_limits<double>::infinity();

        for (const auto &candidate : roots)
        {
            const double z = scalarValue(candidate);
            if (!std::isfinite(z) || z < scalarValue(mix.B))
                continue;
            const auto phase = fugacityAtRoot(pressure, composition, parameters, mix, candidate);
            double g = 0.0;
            for (int i = 0; i < numComponents; ++i)
            {
                const double xi = scalarValue(composition[static_cast<std::size_t>(i)]);
                if (xi > 0.0)
                    g += xi * std::log(std::max(scalarValue(phase.fugacityCoefficient[static_cast<std::size_t>(i)]), 1.0e-300));
            }
            if (!found || g < bestGibbs)
            {
                bestGibbs = g;
                best = candidate;
                found = true;
            }
        }
        if (!found)
            throw std::runtime_error("Cubic EOS Gibbs selection found no admissible root.");
        return best;
    }

    template <class Scalar>
    static void normalizeComposition(Composition<Scalar> &composition, double floor)
    {
        Scalar sum = 0.0;
        for (auto &entry : composition)
        {
            if (scalarValue(entry) < floor)
                entry = entry + (floor - scalarValue(entry));
            sum += entry;
        }
        if (!(scalarValue(sum) > 0.0))
            throw std::runtime_error("Composition normalization requires positive sum.");
        for (auto &entry : composition)
            entry /= sum;
    }

    static constexpr double cpaGasConstant_ = units::gasConstant;

    struct CpaTemperatureCache_
    {
        double temperature{0.0};
        Composition<double> ai{};
        Composition<double> bi{};
        InteractionMatrix<double> aij{};
        InteractionMatrix<double> associationTemperatureFactor{};
        InteractionMatrix<double> associationVolume{};
    };

    template <class Scalar>
    struct CpaMixingState_
    {
        Composition<Scalar> ai{};
        Composition<Scalar> bi{};
        Composition<Scalar> Si{};
        InteractionMatrix<Scalar> aij{};
        InteractionMatrix<Scalar> associationTemperatureFactor{};
        InteractionMatrix<Scalar> associationVolume{};
        Scalar a{0.0};
        Scalar b{0.0};
    };

    template <class Scalar>
    struct CpaAssociationState_
    {
        Composition<Scalar> donorX{};
        Composition<Scalar> acceptorX{};
        Scalar g{1.0};
        Scalar gPrimeOverG{0.0};
        Scalar bondedSiteSum{0.0};
    };

    template <class Scalar, std::size_t Size>
    [[nodiscard]] static bool solveCpaDenseSystem_(
        std::array<std::array<Scalar, Size>, Size> &matrix,
        std::array<Scalar, Size> &rightHandSide,
        int activeSize,
        std::array<Scalar, Size> &solution)
    {
        for (int column = 0; column < activeSize; ++column)
        {
            int pivot = column;
            double pivotMagnitude =
                std::abs(scalarValue(matrix[static_cast<std::size_t>(column)]
                                           [static_cast<std::size_t>(column)]));
            for (int row = column + 1; row < activeSize; ++row)
            {
                const double magnitude =
                    std::abs(scalarValue(matrix[static_cast<std::size_t>(row)]
                                               [static_cast<std::size_t>(column)]));
                if (magnitude > pivotMagnitude)
                {
                    pivot = row;
                    pivotMagnitude = magnitude;
                }
            }
            if (!std::isfinite(pivotMagnitude) || pivotMagnitude <= 1.0e-18)
                return false;

            if (pivot != column)
            {
                std::swap(matrix[static_cast<std::size_t>(pivot)],
                          matrix[static_cast<std::size_t>(column)]);
                std::swap(rightHandSide[static_cast<std::size_t>(pivot)],
                          rightHandSide[static_cast<std::size_t>(column)]);
            }

            const Scalar diagonal =
                matrix[static_cast<std::size_t>(column)]
                      [static_cast<std::size_t>(column)];
            for (int row = column + 1; row < activeSize; ++row)
            {
                const std::size_t rr = static_cast<std::size_t>(row);
                const std::size_t cc = static_cast<std::size_t>(column);
                const Scalar factor = matrix[rr][cc] / diagonal;
                matrix[rr][cc] = Scalar(0.0);
                for (int j = column + 1; j < activeSize; ++j)
                {
                    const std::size_t jj = static_cast<std::size_t>(j);
                    matrix[rr][jj] -= factor * matrix[cc][jj];
                }
                rightHandSide[rr] -= factor * rightHandSide[cc];
            }
        }

        for (int row = activeSize - 1; row >= 0; --row)
        {
            const std::size_t rr = static_cast<std::size_t>(row);
            Scalar value = rightHandSide[rr];
            for (int column = row + 1; column < activeSize; ++column)
            {
                const std::size_t cc = static_cast<std::size_t>(column);
                value -= matrix[rr][cc] * solution[cc];
            }
            const Scalar diagonal = matrix[rr][rr];
            const double diagonalValue = std::abs(scalarValue(diagonal));
            if (!std::isfinite(diagonalValue) || diagonalValue <= 1.0e-18)
                return false;
            solution[rr] = value / diagonal;
            if (!std::isfinite(scalarValue(solution[rr])))
                return false;
        }
        return true;
    }

    [[nodiscard]] bool cpaWaterOnlyAnalyticFastPathAvailable_() const noexcept
    {
        if (!cubicPlusAssociation_)
            return false;
        const auto &options = *cubicPlusAssociation_;
        int associatingComponent = -1;
        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            const int donor = options.donorSites[idx];
            const int acceptor = options.acceptorSites[idx];
            if (donor == 0 && acceptor == 0)
                continue;
            if (associatingComponent >= 0 || donor <= 0 || donor != acceptor)
                return false;
            associatingComponent = i;
        }
        return associatingComponent >= 0;
    }

    [[nodiscard]] int cpaSingleSymmetricAssociatingComponent_() const noexcept
    {
        if (!cpaWaterOnlyAnalyticFastPathAvailable_())
            return -1;
        const auto &options = *cubicPlusAssociation_;
        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            if (options.donorSites[idx] > 0)
                return i;
        }
        return -1;
    }

    [[nodiscard]] int cpaSingleDonorComponent_() const noexcept
    {
        if (!cubicPlusAssociation_)
            return -1;
        const auto &options = *cubicPlusAssociation_;
        int donorComponent = -1;
        for (int i = 0; i < numComponents; ++i)
        {
            if (options.donorSites[static_cast<std::size_t>(i)] <= 0)
                continue;
            if (donorComponent >= 0)
                return -1;
            donorComponent = i;
        }
        return donorComponent;
    }

    [[nodiscard]] CpaTemperatureCache_ buildCpaTemperatureCache_(double temperature) const
    {
        const auto &options = cpaOptions_();
        CpaTemperatureCache_ cache;
        cache.temperature = temperature;
        Composition<double> sqrtAi{};

        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            const double tr = temperature / mixture_.criticalTemperature(i);
            const double alphaRoot = 1.0 + options.c1[idx] * (1.0 - std::sqrt(tr));
            cache.ai[idx] = options.a0[idx] * alphaRoot * alphaRoot;
            cache.bi[idx] = options.b[idx];
            sqrtAi[idx] = std::sqrt(cache.ai[idx]);
        }

        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t ii = static_cast<std::size_t>(i);
            for (int j = 0; j < numComponents; ++j)
            {
                const std::size_t jj = static_cast<std::size_t>(j);
                cache.aij[ii][jj] = sqrtAi[ii] * sqrtAi[jj] *
                    (1.0 - binaryInteractionCoefficient(i, j, temperature));

                if (options.donorSites[ii] <= 0 || options.acceptorSites[jj] <= 0)
                    continue;

                double epsilon = options.crossAssociationEnergy[ii][jj];
                double beta = options.crossAssociationVolume[ii][jj];
                if (!(epsilon > 0.0) || !(beta > 0.0))
                {
                    if (options.associationEnergy[ii] > 0.0 &&
                        options.associationEnergy[jj] > 0.0 &&
                        options.associationVolume[ii] > 0.0 &&
                        options.associationVolume[jj] > 0.0)
                    {
                        epsilon = 0.5 *
                            (options.associationEnergy[ii] + options.associationEnergy[jj]);
                        beta = std::sqrt(
                            options.associationVolume[ii] * options.associationVolume[jj]);
                    }
                    else
                    {
                        epsilon = 0.0;
                        beta = 0.0;
                    }
                }
                if (epsilon > 0.0 && beta > 0.0)
                {
                    cache.associationTemperatureFactor[ii][jj] =
                        std::exp(epsilon / (cpaGasConstant_ * temperature)) - 1.0;
                    cache.associationVolume[ii][jj] = beta;
                }
            }
        }
        return cache;
    }

    [[nodiscard]] const CubicPlusAssociationOptions &cpaOptions_() const
    {
        if (!cubicPlusAssociation_)
            throw std::logic_error("CPA options are not configured.");
        return *cubicPlusAssociation_;
    }

    template <class Scalar>
    [[nodiscard]] CpaMixingState_<Scalar> cpaMixing_(
        const Scalar &temperature,
        const Composition<Scalar> &composition) const
    {
        const bool profiling = thermodynamicProfilerEnabled_;
        const auto profileStart = profiling ? std::chrono::steady_clock::now()
                                            : std::chrono::steady_clock::time_point{};
        if (profiling)
            ++thermodynamicProfile_.cpaMixingCalls;

        const auto &options = cpaOptions_();
        CpaMixingState_<Scalar> result;
        const double tValue = scalarValue(temperature);
        const bool cacheHit = cpaTemperatureCache_.has_value() &&
            std::abs(tValue - cpaTemperatureCache_->temperature) <=
                1.0e-12 * std::max(1.0, std::abs(tValue));

        if (cacheHit)
        {
            const auto &cache = *cpaTemperatureCache_;
            if (profiling)
                ++thermodynamicProfile_.cpaTemperatureCacheHits;
            for (int i = 0; i < numComponents; ++i)
            {
                const std::size_t ii = static_cast<std::size_t>(i);
                result.ai[ii] = Scalar(cache.ai[ii]);
                result.bi[ii] = Scalar(cache.bi[ii]);
                for (int j = 0; j < numComponents; ++j)
                {
                    const std::size_t jj = static_cast<std::size_t>(j);
                    result.aij[ii][jj] = Scalar(cache.aij[ii][jj]);
                    result.associationTemperatureFactor[ii][jj] =
                        Scalar(cache.associationTemperatureFactor[ii][jj]);
                    result.associationVolume[ii][jj] =
                        Scalar(cache.associationVolume[ii][jj]);
                }
            }
        }
        else
        {
            Composition<Scalar> sqrtAi{};
            for (int i = 0; i < numComponents; ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(i);
                const Scalar tr = temperature / mixture_.criticalTemperature(i);
                const Scalar alphaRoot = Scalar(1.0) + options.c1[idx] *
                    (Scalar(1.0) - sqrt(tr));
                result.ai[idx] = options.a0[idx] * alphaRoot * alphaRoot;
                result.bi[idx] = options.b[idx];
                sqrtAi[idx] = sqrt(result.ai[idx]);
            }

            for (int i = 0; i < numComponents; ++i)
            {
                const std::size_t ii = static_cast<std::size_t>(i);
                for (int j = 0; j < numComponents; ++j)
                {
                    const std::size_t jj = static_cast<std::size_t>(j);
                    const double kij = binaryInteractionCoefficient(i, j, tValue);
                    result.aij[ii][jj] = sqrtAi[ii] * sqrtAi[jj] * (1.0 - kij);

                    if (options.donorSites[ii] > 0 && options.acceptorSites[jj] > 0)
                    {
                        double epsilon = options.crossAssociationEnergy[ii][jj];
                        double beta = options.crossAssociationVolume[ii][jj];
                        if (!(epsilon > 0.0) || !(beta > 0.0))
                        {
                            if (options.associationEnergy[ii] > 0.0 &&
                                options.associationEnergy[jj] > 0.0 &&
                                options.associationVolume[ii] > 0.0 &&
                                options.associationVolume[jj] > 0.0)
                            {
                                epsilon = 0.5 *
                                    (options.associationEnergy[ii] + options.associationEnergy[jj]);
                                beta = std::sqrt(
                                    options.associationVolume[ii] *
                                    options.associationVolume[jj]);
                            }
                            else
                            {
                                epsilon = 0.0;
                                beta = 0.0;
                            }
                        }

                        if (epsilon > 0.0 && beta > 0.0)
                        {
                            result.associationTemperatureFactor[ii][jj] =
                                exp(epsilon / (cpaGasConstant_ * temperature)) - 1.0;
                            result.associationVolume[ii][jj] = beta;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t ii = static_cast<std::size_t>(i);
            for (int j = 0; j < numComponents; ++j)
            {
                const std::size_t jj = static_cast<std::size_t>(j);
                result.Si[ii] += composition[jj] * result.aij[ii][jj];
            }
            result.a += composition[ii] * result.Si[ii];
            result.b += composition[ii] * result.bi[ii];
        }

        if (profiling)
        {
            thermodynamicProfile_.cpaMixingSeconds +=
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - profileStart).count();
        }
        return result;
    }

    template <class Scalar>
    [[nodiscard]] std::pair<Scalar, Scalar> cpaRadialDistribution_(
        const Scalar &eta) const
    {
        const auto radial = cpaOptions_().radialDistribution;
        if (radial == CpaRadialDistribution::Simplified)
        {
            const Scalar denominator = Scalar(1.0) - 1.9 * eta;
            const Scalar g = Scalar(1.0) / denominator;
            return {g, 1.9 * g}; // g'/g
        }

        const Scalar oneMinusEta = Scalar(1.0) - eta;
        const Scalar g = (Scalar(1.0) - 0.5 * eta) /
            (oneMinusEta * oneMinusEta * oneMinusEta);
        const Scalar gPrime = (Scalar(2.5) - eta) /
            (oneMinusEta * oneMinusEta * oneMinusEta * oneMinusEta);
        return {g, gPrime / g};
    }

    template <class Scalar>
    [[nodiscard]] Scalar cpaAssociationDelta_(
        int donorComponent,
        int acceptorComponent,
        const CpaMixingState_<Scalar> &mix,
        const Scalar &g) const
    {
        const auto i = static_cast<std::size_t>(donorComponent);
        const auto j = static_cast<std::size_t>(acceptorComponent);
        const Scalar temperatureFactor = mix.associationTemperatureFactor[i][j];
        const Scalar beta = mix.associationVolume[i][j];
        if (!(scalarValue(temperatureFactor) > 0.0) ||
            !(scalarValue(beta) > 0.0))
        {
            return Scalar(0.0);
        }

        const Scalar bij = 0.5 * (mix.bi[i] + mix.bi[j]);
        // 保持原来的 g*(exp-1)*bij*beta 乘法顺序，使 cache 只消除重复的
        // 超越函数和温度 mixing 工作，不改变原公式的数值顺序。
        return g * temperatureFactor * bij * beta;
    }

    template <class Scalar>
    [[nodiscard]] CpaAssociationState_<Scalar> cpaAssociationState_(
        const Scalar &rho,
        const Composition<Scalar> &composition,
        const CpaMixingState_<Scalar> &mix) const
    {
        const bool profiling = thermodynamicProfilerEnabled_;
        const auto profileStart = profiling ? std::chrono::steady_clock::now()
                                            : std::chrono::steady_clock::time_point{};
        if (profiling)
            ++thermodynamicProfile_.cpaAssociationCalls;

        const auto &options = cpaOptions_();
        CpaAssociationState_<Scalar> state;
        for (auto &entry : state.donorX) entry = Scalar(1.0);
        for (auto &entry : state.acceptorX) entry = Scalar(1.0);

        const Scalar eta = mix.b * rho / 4.0;
        auto radial = cpaRadialDistribution_(eta);
        state.g = radial.first;
        state.gPrimeOverG = radial.second;

        // 常见 4C-water 快速路径：只有一个缔合组分，且 donor/acceptor 位点数
        // 相同。由对称性 X_D=X_A=X，质量作用方程可化为
        // A X^2 + X - 1 = 0，其中 A=rho*x*nSites*Delta。
        const int symmetricComponent = cpaSingleSymmetricAssociatingComponent_();
        if (symmetricComponent >= 0)
        {
            const std::size_t k = static_cast<std::size_t>(symmetricComponent);
            const int sites = options.donorSites[k];
            const Scalar delta = cpaAssociationDelta_(
                symmetricComponent, symmetricComponent, mix, state.g);
            const Scalar A = rho * composition[k] * static_cast<double>(sites) * delta;
            const Scalar X = Scalar(2.0) /
                (Scalar(1.0) + sqrt(Scalar(1.0) + Scalar(4.0) * A));
            state.donorX[k] = X;
            state.acceptorX[k] = X;
            if (profiling)
                ++thermodynamicProfile_.cpaAssociationAnalyticCalls;
        }
        else if (const int donorComponent = cpaSingleDonorComponent_();
                 donorComponent >= 0)
        {
            // A single donor-bearing component (for example 4C water in the
            // B2 H2O-CO2 scheme) reduces the coupled site equations exactly to
            // one scalar equation.  For donor fraction X_D,
            //
            //   X_A,j = 1 / (1 + B_j X_D)
            //   X_D (1 + sum_j A_j / (1 + B_j X_D)) - 1 = 0,
            //
            // where A_j=rho*x_j*n_A,j*Delta_Dj and
            // B_j=rho*x_D*n_D*Delta_Dj.  A safeguarded Newton iteration on
            // this monotone equation replaces hundreds of damped updates while
            // retaining Scalar arithmetic so automatic derivatives propagate.
            if (profiling)
                ++thermodynamicProfile_.cpaAssociationIterativeCalls;

            const std::size_t donor = static_cast<std::size_t>(donorComponent);
            Composition<Scalar> aCoefficient{};
            Composition<Scalar> bCoefficient{};
            Scalar aSum = 0.0;
            for (int j = 0; j < numComponents; ++j)
            {
                const std::size_t acceptor = static_cast<std::size_t>(j);
                if (options.acceptorSites[acceptor] <= 0)
                    continue;
                const Scalar delta = cpaAssociationDelta_(
                    donorComponent, j, mix, state.g);
                aCoefficient[acceptor] = rho * composition[acceptor] *
                    options.acceptorSites[acceptor] * delta;
                bCoefficient[acceptor] = rho * composition[donor] *
                    options.donorSites[donor] * delta;
                aSum += aCoefficient[acceptor];
            }

            Scalar donorX = Scalar(1.0) / (Scalar(1.0) + aSum);
            double lower = 0.0;
            double upper = 1.0;
            bool converged = false;
            const int maximumIterations = std::min(
                options.maximumAssociationIterations, 64);
            for (int iteration = 0; iteration < maximumIterations; ++iteration)
            {
                if (profiling)
                    ++thermodynamicProfile_.cpaAssociationIterations;

                Scalar sum = 0.0;
                Scalar slope = 0.0;
                for (int j = 0; j < numComponents; ++j)
                {
                    const std::size_t acceptor = static_cast<std::size_t>(j);
                    if (options.acceptorSites[acceptor] <= 0)
                        continue;
                    const Scalar denominator =
                        Scalar(1.0) + bCoefficient[acceptor] * donorX;
                    sum += aCoefficient[acceptor] / denominator;
                    slope -= aCoefficient[acceptor] * bCoefficient[acceptor] /
                        (denominator * denominator);
                }

                const Scalar residual = donorX * (Scalar(1.0) + sum) - 1.0;
                const Scalar derivative =
                    Scalar(1.0) + sum + donorX * slope;
                const double residualValue = scalarValue(residual);
                const double donorValue = scalarValue(donorX);
                if (residualValue > 0.0)
                    upper = donorValue;
                else
                    lower = donorValue;

                if (std::abs(residualValue) <= options.associationTolerance)
                {
                    converged = true;
                    break;
                }

                const double derivativeValue = scalarValue(derivative);
                Scalar candidate = donorX;
                if (std::isfinite(derivativeValue) &&
                    std::abs(derivativeValue) > 1.0e-14)
                {
                    candidate = donorX - residual / derivative;
                }
                const double candidateValue = scalarValue(candidate);
                if (!std::isfinite(candidateValue) ||
                    !(candidateValue > lower && candidateValue < upper))
                {
                    const double midpoint = 0.5 * (lower + upper);
                    candidate = donorX + (midpoint - donorValue);
                }
                donorX = candidate;
            }
            if (!converged)
                throw std::runtime_error(
                    "CPA reduced association-site iteration did not converge.");

            state.donorX[donor] = donorX;
            for (int j = 0; j < numComponents; ++j)
            {
                const std::size_t acceptor = static_cast<std::size_t>(j);
                if (options.acceptorSites[acceptor] > 0)
                {
                    state.acceptorX[acceptor] = Scalar(1.0) /
                        (Scalar(1.0) + bCoefficient[acceptor] * donorX);
                }
            }
        }
        else
        {
            if (profiling)
                ++thermodynamicProfile_.cpaAssociationIterativeCalls;

            // General multi-donor/multi-acceptor path.  Its mass-action
            // residuals are small (at most 2*numComponents), so a dense Newton
            // step is substantially cheaper than hundreds of damped fixed-point
            // sweeps.  Scalar-valued elimination preserves automatic
            // derivatives.  A bounded line search enforces 0<X<=1, and the
            // historical fixed-point iteration remains as a fail-safe.
            constexpr std::size_t maximumUnknowns =
                static_cast<std::size_t>(2 * numComponents);
            std::array<int, numComponents> donorUnknown{};
            std::array<int, numComponents> acceptorUnknown{};
            donorUnknown.fill(-1);
            acceptorUnknown.fill(-1);
            int unknownCount = 0;
            for (int i = 0; i < numComponents; ++i)
            {
                if (options.donorSites[static_cast<std::size_t>(i)] > 0)
                    donorUnknown[static_cast<std::size_t>(i)] = unknownCount++;
            }
            for (int i = 0; i < numComponents; ++i)
            {
                if (options.acceptorSites[static_cast<std::size_t>(i)] > 0)
                    acceptorUnknown[static_cast<std::size_t>(i)] = unknownCount++;
            }

            InteractionMatrix<Scalar> donorCoupling{};
            InteractionMatrix<Scalar> acceptorCoupling{};
            for (int i = 0; i < numComponents; ++i)
            {
                const std::size_t donor = static_cast<std::size_t>(i);
                if (options.donorSites[donor] <= 0)
                    continue;
                for (int j = 0; j < numComponents; ++j)
                {
                    const std::size_t acceptor = static_cast<std::size_t>(j);
                    if (options.acceptorSites[acceptor] <= 0)
                        continue;
                    const Scalar delta = cpaAssociationDelta_(i, j, mix, state.g);
                    donorCoupling[donor][acceptor] = rho * composition[acceptor] *
                        options.acceptorSites[acceptor] * delta;
                    acceptorCoupling[acceptor][donor] = rho * composition[donor] *
                        options.donorSites[donor] * delta;
                }
            }

            // One Jacobi mass-action update from X=1 gives a physical and much
            // tighter initial point than the unassociated state.
            for (int i = 0; i < numComponents; ++i)
            {
                const std::size_t component = static_cast<std::size_t>(i);
                if (options.donorSites[component] > 0)
                {
                    Scalar sum = 0.0;
                    for (int j = 0; j < numComponents; ++j)
                        sum += donorCoupling[component][static_cast<std::size_t>(j)];
                    state.donorX[component] = Scalar(1.0) / (Scalar(1.0) + sum);
                }
                if (options.acceptorSites[component] > 0)
                {
                    Scalar sum = 0.0;
                    for (int j = 0; j < numComponents; ++j)
                        sum += acceptorCoupling[component][static_cast<std::size_t>(j)];
                    state.acceptorX[component] = Scalar(1.0) / (Scalar(1.0) + sum);
                }
            }

            const auto residualNorm = [&](const Composition<Scalar> &donorX,
                                          const Composition<Scalar> &acceptorX) {
                double norm = 0.0;
                for (int i = 0; i < numComponents; ++i)
                {
                    const std::size_t component = static_cast<std::size_t>(i);
                    if (options.donorSites[component] > 0)
                    {
                        Scalar denominator = 1.0;
                        for (int j = 0; j < numComponents; ++j)
                        {
                            const std::size_t acceptor = static_cast<std::size_t>(j);
                            denominator += donorCoupling[component][acceptor] *
                                acceptorX[acceptor];
                        }
                        norm = std::max(norm, std::abs(scalarValue(
                            donorX[component] * denominator - 1.0)));
                    }
                    if (options.acceptorSites[component] > 0)
                    {
                        Scalar denominator = 1.0;
                        for (int j = 0; j < numComponents; ++j)
                        {
                            const std::size_t donor = static_cast<std::size_t>(j);
                            denominator += acceptorCoupling[component][donor] *
                                donorX[donor];
                        }
                        norm = std::max(norm, std::abs(scalarValue(
                            acceptorX[component] * denominator - 1.0)));
                    }
                }
                return norm;
            };

            bool newtonConverged = false;
            const int maximumNewtonIterations = std::min(
                options.maximumAssociationIterations, 64);
            for (int iteration = 0;
                 iteration < maximumNewtonIterations;
                 ++iteration)
            {
                if (profiling)
                    ++thermodynamicProfile_.cpaAssociationIterations;

                std::array<std::array<Scalar, maximumUnknowns>, maximumUnknowns>
                    jacobian{};
                std::array<Scalar, maximumUnknowns> rightHandSide{};
                std::array<Scalar, maximumUnknowns> increment{};
                double norm = 0.0;

                for (int i = 0; i < numComponents; ++i)
                {
                    const std::size_t component = static_cast<std::size_t>(i);
                    const int row = donorUnknown[component];
                    if (row < 0)
                        continue;
                    Scalar denominator = 1.0;
                    for (int j = 0; j < numComponents; ++j)
                    {
                        const std::size_t acceptor = static_cast<std::size_t>(j);
                        denominator += donorCoupling[component][acceptor] *
                            state.acceptorX[acceptor];
                        const int column = acceptorUnknown[acceptor];
                        if (column >= 0)
                        {
                            jacobian[static_cast<std::size_t>(row)]
                                    [static_cast<std::size_t>(column)] =
                                state.donorX[component] *
                                donorCoupling[component][acceptor];
                        }
                    }
                    const Scalar residual =
                        state.donorX[component] * denominator - 1.0;
                    rightHandSide[static_cast<std::size_t>(row)] = -residual;
                    jacobian[static_cast<std::size_t>(row)]
                            [static_cast<std::size_t>(row)] = denominator;
                    norm = std::max(norm, std::abs(scalarValue(residual)));
                }

                for (int i = 0; i < numComponents; ++i)
                {
                    const std::size_t component = static_cast<std::size_t>(i);
                    const int row = acceptorUnknown[component];
                    if (row < 0)
                        continue;
                    Scalar denominator = 1.0;
                    for (int j = 0; j < numComponents; ++j)
                    {
                        const std::size_t donor = static_cast<std::size_t>(j);
                        denominator += acceptorCoupling[component][donor] *
                            state.donorX[donor];
                        const int column = donorUnknown[donor];
                        if (column >= 0)
                        {
                            jacobian[static_cast<std::size_t>(row)]
                                    [static_cast<std::size_t>(column)] =
                                state.acceptorX[component] *
                                acceptorCoupling[component][donor];
                        }
                    }
                    const Scalar residual =
                        state.acceptorX[component] * denominator - 1.0;
                    rightHandSide[static_cast<std::size_t>(row)] = -residual;
                    jacobian[static_cast<std::size_t>(row)]
                            [static_cast<std::size_t>(row)] = denominator;
                    norm = std::max(norm, std::abs(scalarValue(residual)));
                }

                if (norm <= options.associationTolerance)
                {
                    newtonConverged = true;
                    break;
                }
                if (!solveCpaDenseSystem_(
                        jacobian, rightHandSide, unknownCount, increment))
                {
                    break;
                }

                bool accepted = false;
                double damping = 1.0;
                for (int lineSearch = 0; lineSearch < 20; ++lineSearch)
                {
                    Composition<Scalar> candidateDonor = state.donorX;
                    Composition<Scalar> candidateAcceptor = state.acceptorX;
                    bool physical = true;
                    for (int i = 0; i < numComponents; ++i)
                    {
                        const std::size_t component = static_cast<std::size_t>(i);
                        if (donorUnknown[component] >= 0)
                        {
                            candidateDonor[component] += damping * increment[
                                static_cast<std::size_t>(donorUnknown[component])];
                            const double value = scalarValue(candidateDonor[component]);
                            physical = physical && std::isfinite(value) &&
                                value > 0.0 && value <= 1.0;
                        }
                        if (acceptorUnknown[component] >= 0)
                        {
                            candidateAcceptor[component] += damping * increment[
                                static_cast<std::size_t>(acceptorUnknown[component])];
                            const double value = scalarValue(candidateAcceptor[component]);
                            physical = physical && std::isfinite(value) &&
                                value > 0.0 && value <= 1.0;
                        }
                    }

                    if (physical)
                    {
                        const double candidateNorm =
                            residualNorm(candidateDonor, candidateAcceptor);
                        if (std::isfinite(candidateNorm) && candidateNorm < norm)
                        {
                            state.donorX = std::move(candidateDonor);
                            state.acceptorX = std::move(candidateAcceptor);
                            accepted = true;
                            if (candidateNorm <= options.associationTolerance)
                                newtonConverged = true;
                            break;
                        }
                    }
                    damping *= 0.5;
                }
                if (newtonConverged || !accepted)
                    break;
            }

            if (!newtonConverged)
            {
                for (auto &entry : state.donorX) entry = Scalar(1.0);
                for (auto &entry : state.acceptorX) entry = Scalar(1.0);
                for (int iteration = 0;
                     iteration < options.maximumAssociationIterations;
                     ++iteration)
                {
                    if (profiling)
                        ++thermodynamicProfile_.cpaAssociationIterations;
                    Composition<Scalar> newDonor = state.donorX;
                    Composition<Scalar> newAcceptor = state.acceptorX;
                    double maxChange = 0.0;

                    for (int i = 0; i < numComponents; ++i)
                    {
                        const std::size_t ii = static_cast<std::size_t>(i);
                        if (options.donorSites[ii] > 0)
                        {
                            Scalar sum = 0.0;
                            for (int j = 0; j < numComponents; ++j)
                            {
                                const std::size_t jj = static_cast<std::size_t>(j);
                                sum += donorCoupling[ii][jj] * state.acceptorX[jj];
                            }
                            const Scalar target = Scalar(1.0) /
                                (Scalar(1.0) + sum);
                            newDonor[ii] = state.donorX[ii] +
                                options.associationDamping *
                                (target - state.donorX[ii]);
                            maxChange = std::max(maxChange,
                                std::abs(scalarValue(newDonor[ii] - state.donorX[ii])));
                        }
                        if (options.acceptorSites[ii] > 0)
                        {
                            Scalar sum = 0.0;
                            for (int j = 0; j < numComponents; ++j)
                            {
                                const std::size_t jj = static_cast<std::size_t>(j);
                                sum += acceptorCoupling[ii][jj] * state.donorX[jj];
                            }
                            const Scalar target = Scalar(1.0) /
                                (Scalar(1.0) + sum);
                            newAcceptor[ii] = state.acceptorX[ii] +
                                options.associationDamping *
                                (target - state.acceptorX[ii]);
                            maxChange = std::max(maxChange,
                                std::abs(scalarValue(newAcceptor[ii] - state.acceptorX[ii])));
                        }
                    }

                    state.donorX = std::move(newDonor);
                    state.acceptorX = std::move(newAcceptor);
                    if (maxChange <= options.associationTolerance)
                        break;

                    if (iteration + 1 == options.maximumAssociationIterations)
                    {
                        throw std::runtime_error(
                            "CPA association-site iteration did not converge.");
                    }
                }
            }
        }

        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            state.bondedSiteSum += composition[idx] *
                (options.donorSites[idx] * (Scalar(1.0) - state.donorX[idx]) +
                 options.acceptorSites[idx] * (Scalar(1.0) - state.acceptorX[idx]));
        }

        if (profiling)
        {
            thermodynamicProfile_.cpaAssociationSeconds +=
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - profileStart).count();
        }
        return state;
    }

    template <class Scalar>
    [[nodiscard]] Scalar cpaPressureAtDensityWithMix_(
        const Scalar &rho,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        const CpaMixingState_<Scalar> &mix) const
    {
        if (!(scalarValue(Scalar(1.0) - mix.b * rho) > 0.0))
            throw std::runtime_error("CPA density exceeds the SRK co-volume limit.");

        const Scalar q = mix.b * rho;
        Scalar attractionDenominator = Scalar(1.0) + q;
        if (cpaOptions_().physicalTerm == CpaCubicPhysicalTerm::PengRobinson)
        {
            constexpr double sqrtTwo = 1.4142135623730950488;
            attractionDenominator =
                (Scalar(1.0) + (1.0 + sqrtTwo) * q) *
                (Scalar(1.0) + (1.0 - sqrtTwo) * q);
        }
        const Scalar cubicPressure = cpaGasConstant_ * temperature * rho /
                (Scalar(1.0) - q)
            - mix.a * rho * rho / attractionDenominator;
        const auto assoc = cpaAssociationState_(rho, composition, mix);
        const Scalar factor = Scalar(1.0) +
            assoc.gPrimeOverG * mix.b * rho / 4.0;
        const Scalar associationPressure = -0.5 * cpaGasConstant_ * temperature *
            rho * factor * assoc.bondedSiteSum;
        return cubicPressure + associationPressure;
    }

    template <class Scalar>
    [[nodiscard]] Scalar cpaPressureAtDensity_(
        const Scalar &rho,
        const Scalar &temperature,
        const Composition<Scalar> &composition) const
    {
        const auto mix = cpaMixing_(temperature, composition);
        return cpaPressureAtDensityWithMix_(rho, temperature, composition, mix);
    }

    [[nodiscard]] double cpaDensityRootValueWithMix_(
        double pressure,
        double temperature,
        const Composition<double> &composition,
        bool liquid,
        const CpaMixingState_<double> &mix) const
    {
        const bool profiling = thermodynamicProfilerEnabled_;
        const auto profileStart = profiling ? std::chrono::steady_clock::now()
                                            : std::chrono::steady_clock::time_point{};
        if (profiling)
            ++thermodynamicProfile_.cpaDensityRootCalls;

        const double b = scalarValue(mix.b);
        if (!(b > 0.0))
            throw std::runtime_error("CPA mixture co-volume must be positive.");
        const double rhoMax = 0.999999 / b;
        auto residual = [&](double rho) {
            if (profiling)
                ++thermodynamicProfile_.cpaDensityResidualEvaluations;
            return scalarValue(cpaPressureAtDensityWithMix_(
                rho, temperature, composition, mix)) - pressure;
        };

        double lo = 0.0;
        double hi = 0.0;
        double flo = -pressure;
        double fhi = 0.0;

        if (!liquid)
        {
            hi = std::min(rhoMax, std::max(1.0, 0.25 * pressure /
                (cpaGasConstant_ * temperature)));
            fhi = residual(hi);
            while (fhi < 0.0 && hi < 0.999 * rhoMax)
            {
                lo = hi;
                flo = fhi;
                hi = std::min(0.999999 * rhoMax, hi * 1.6 + 1.0);
                fhi = residual(hi);
            }
        }
        else
        {
            hi = rhoMax;
            fhi = residual(hi);
            lo = 0.82 * hi;
            flo = residual(lo);
            for (int i = 0; i < 160 && flo > 0.0; ++i)
            {
                hi = lo;
                fhi = flo;
                lo *= 0.82;
                flo = residual(lo);
            }
        }

        if (!(flo <= 0.0 && fhi >= 0.0))
            throw std::runtime_error("CPA could not bracket the requested density root.");

        for (int iteration = 0; iteration < 100; ++iteration)
        {
            if (profiling)
                ++thermodynamicProfile_.cpaDensityBisectionIterations;
            const double mid = 0.5 * (lo + hi);
            const double fm = residual(mid);
            if (fm > 0.0)
            {
                hi = mid;
                fhi = fm;
            }
            else
            {
                lo = mid;
                flo = fm;
            }
            if (std::abs(fm) <= 1.0e-10 * std::max(pressure, 1.0) ||
                (hi - lo) <= 1.0e-12 * std::max(mid, 1.0))
                break;
        }
        const double root = 0.5 * (lo + hi);
        if (profiling)
        {
            thermodynamicProfile_.cpaDensityRootSeconds +=
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - profileStart).count();
        }
        return root;
    }

    [[nodiscard]] double cpaDensityRootValue_(
        double pressure,
        double temperature,
        const Composition<double> &composition,
        bool liquid) const
    {
        const auto mix = cpaMixing_(temperature, composition);
        return cpaDensityRootValueWithMix_(
            pressure, temperature, composition, liquid, mix);
    }

    template <class Scalar>
    [[nodiscard]] Scalar cpaDensityRoot_(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        bool liquid) const
    {
        Composition<double> scalarComposition{};
        for (int i = 0; i < numComponents; ++i)
            scalarComposition[static_cast<std::size_t>(i)] =
                scalarValue(composition[static_cast<std::size_t>(i)]);

        const double pValue = scalarValue(pressure);
        const double tValue = scalarValue(temperature);
        const auto scalarMix = cpaMixing_(tValue, scalarComposition);
        const double root = cpaDensityRootValueWithMix_(
            pValue, tValue, scalarComposition, liquid, scalarMix);

        Scalar rho = pressure;
        rho *= 0.0;
        rho += root;

        // Inject first derivatives through the implicit pressure equation.
        // One correction is sufficient for first-order AD at an already converged
        // scalar root; a second correction also removes the scalar root tolerance.
        const auto applyCorrections = [&](const auto &differentiableMix) {
            for (int correction = 0; correction < 2; ++correction)
            {
                const double rhoValue = scalarValue(rho);
                const double h = std::max(1.0e-4, std::abs(rhoValue) * 2.0e-6);
                const double lower = std::max(1.0e-12, rhoValue - h);
                const double upper = rhoValue + h;
                const double pLower = scalarValue(cpaPressureAtDensityWithMix_(
                    lower, tValue, scalarComposition, scalarMix));
                const double pUpper = scalarValue(cpaPressureAtDensityWithMix_(
                    upper, tValue, scalarComposition, scalarMix));
                const double derivative = (pUpper - pLower) / (upper - lower);
                if (!std::isfinite(derivative) || std::abs(derivative) < 1.0e-12)
                    throw std::runtime_error("CPA density root has a singular pressure derivative.");
                const Scalar residual = cpaPressureAtDensityWithMix_(
                    rho, temperature, composition, differentiableMix) - pressure;
                rho -= residual / derivative;
            }
        };

        if constexpr (std::is_same_v<Scalar, double>)
        {
            // For scalar property/flash evaluations the already-built scalar
            // mixing state is exactly the differentiable state as well.
            applyCorrections(scalarMix);
        }
        else
        {
            const auto differentiableMix = cpaMixing_(temperature, composition);
            applyCorrections(differentiableMix);
        }
        return rho;
    }

    template <class Scalar>
    [[nodiscard]] PhaseResult<Scalar> cpaAtDensity_(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        const Scalar &rho) const
    {
        const auto mix = cpaMixing_(temperature, composition);
        const auto assoc = cpaAssociationState_(rho, composition, mix);
        const Scalar q = mix.b * rho;
        if (!(scalarValue(Scalar(1.0) - q) > 0.0))
            throw std::runtime_error("CPA fugacity requires b*rho < 1.");

        PhaseResult<Scalar> result;
        result.compressibility = pressure / (rho * cpaGasConstant_ * temperature);
        Scalar cubicHelmholtzIntegral = log(Scalar(1.0) + q);
        Scalar cubicHelmholtzDerivative = Scalar(1.0) /
            (Scalar(1.0) + q);
        if (cpaOptions_().physicalTerm == CpaCubicPhysicalTerm::PengRobinson)
        {
            constexpr double sqrtTwo = 1.4142135623730950488;
            constexpr double deltaOne = 1.0 + sqrtTwo;
            constexpr double deltaTwo = 1.0 - sqrtTwo;
            constexpr double deltaDifference = 2.0 * sqrtTwo;
            const Scalar first = Scalar(1.0) + deltaOne * q;
            const Scalar second = Scalar(1.0) + deltaTwo * q;
            cubicHelmholtzIntegral = log(first / second) / deltaDifference;
            cubicHelmholtzDerivative = Scalar(1.0) / (first * second);
        }
        const Scalar logZ = log(result.compressibility);
        const auto &options = cpaOptions_();

        for (int i = 0; i < numComponents; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            const Scalar bi = mix.bi[idx];
            Scalar muCubic = -log(Scalar(1.0) - q) + rho * bi /
                (Scalar(1.0) - q);
            muCubic -= (
                (2.0 * mix.Si[idx] / mix.b - mix.a * bi / (mix.b * mix.b)) *
                    cubicHelmholtzIntegral
                + (mix.a / mix.b) * rho * bi * cubicHelmholtzDerivative) /
                (cpaGasConstant_ * temperature);

            Scalar muAssociation = 0.0;
            if (options.donorSites[idx] > 0)
                muAssociation += options.donorSites[idx] * log(assoc.donorX[idx]);
            if (options.acceptorSites[idx] > 0)
                muAssociation += options.acceptorSites[idx] * log(assoc.acceptorX[idx]);
            const Scalar dLnGdNi = assoc.gPrimeOverG * bi * rho / 4.0;
            muAssociation -= 0.5 * dLnGdNi * assoc.bondedSiteSum;

            const Scalar lnPhi = muCubic + muAssociation - logZ;
            result.fugacityCoefficient[idx] = exp(lnPhi);
            result.fugacity[idx] = result.fugacityCoefficient[idx] * pressure * composition[idx];
        }
        return result;
    }

    template <class Scalar>
    [[nodiscard]] PhaseResult<Scalar> cpaPhaseResult_(
        const Scalar &pressure,
        const Scalar &temperature,
        const Composition<Scalar> &composition,
        bool liquid,
        bool selectGibbsMinimum) const
    {
        const bool profiling = thermodynamicProfilerEnabled_;
        const auto profileStart = profiling ? std::chrono::steady_clock::now()
                                            : std::chrono::steady_clock::time_point{};
        if (profiling)
            ++thermodynamicProfile_.cpaPhaseResultCalls;
        auto finish = [&](PhaseResult<Scalar> result) {
            if (profiling)
            {
                thermodynamicProfile_.cpaPhaseResultSeconds +=
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - profileStart).count();
            }
            return result;
        };

        const Scalar selectedDensity = cpaDensityRoot_(
            pressure, temperature, composition, liquid);
        auto selected = cpaAtDensity_(
            pressure, temperature, composition, selectedDensity);

        if (!selectGibbsMinimum)
            return finish(std::move(selected));

        try
        {
            const Scalar alternativeDensity = cpaDensityRoot_(
                pressure, temperature, composition, !liquid);
            if (std::abs(scalarValue(alternativeDensity - selectedDensity)) <=
                1.0e-8 * std::max(1.0, std::abs(scalarValue(selectedDensity))))
                return finish(std::move(selected));

            auto alternative = cpaAtDensity_(
                pressure, temperature, composition, alternativeDensity);
            auto reducedGibbs = [&](const PhaseResult<Scalar> &phase) {
                double value = 0.0;
                for (int i = 0; i < numComponents; ++i)
                {
                    const std::size_t idx = static_cast<std::size_t>(i);
                    const double xi = scalarValue(composition[idx]);
                    if (xi > 0.0)
                        value += xi * std::log(std::max(
                            scalarValue(phase.fugacityCoefficient[idx]), 1.0e-300));
                }
                return value;
            };
            if (reducedGibbs(alternative) < reducedGibbs(selected))
                return finish(std::move(alternative));
        }
        catch (const std::runtime_error &)
        {
            // A single physical root is common away from coexistence.
        }
        return finish(std::move(selected));
    }

    std::optional<EquilibriumConstantFunctions> equilibriumConstantFunctions_;
    BinaryInteractionFunction binaryInteractionFunction_{};
    Composition<double> volumeTranslation_{};
    bool volumeTranslationEnabled_{false};
    bool aqueousVolumeTranslationEnabled_{false};
    std::optional<IapwsGarciaAqueousVolume<numComponents>>
        aqueousVolumeModel_{};
    int aqueousDomainWaterComponent_{-1};
    double aqueousDomainMaximumSoluteMoleFraction_{1.0};

    [[nodiscard]] const SoreideWhitsonOptions &soreideWhitsonOptions_() const
    {
        if (!soreideWhitson_)
            throw std::logic_error("Soreide-Whitson options are not configured.");
        return *soreideWhitson_;
    }

    [[nodiscard]] double soreideWhitsonAqueousBip_(
        int i, int j, double temperature) const
    {
        const auto &options = soreideWhitsonOptions_();
        const int w = options.waterComponent;
        if (i != w && j != w)
            return binaryInteractionCoefficient(i, j, temperature);
        if (i == w && j == w)
            return 0.0;
        const int other = (i == w) ? j : i;
        const auto &correlation =
            options.aqueousWaterBip[static_cast<std::size_t>(other)];
        if (correlation)
            return correlation(temperature, options.salinityMolality);
        return binaryInteractionCoefficient(i, j, temperature);
    }

    mutable bool thermodynamicProfilerEnabled_{false};
    mutable ThermodynamicProfile thermodynamicProfile_{};
    std::optional<CpaTemperatureCache_> cpaTemperatureCache_{};

    CubicThermodynamicModel thermodynamicModel_{CubicThermodynamicModel::PengRobinson};
    std::optional<SoreideWhitsonOptions> soreideWhitson_{};
    std::optional<CubicPlusAssociationOptions> cubicPlusAssociation_{};

    double omegaA_{0.4572355};
    double omegaB_{0.0779691};
    CompositionalMixture<Indices> mixture_{};
    int eosType_{1};
    double m1_{2.414213562373095};
    double m2_{-0.414213562373095};
    double minimumComposition_{NaturalNumerics::minimumComposition};
};


} // namespace MPMC
