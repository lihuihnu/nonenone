#pragma once

#include <indices/model_config.hpp>
#include <natural/thermo/cubic_plus_association.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <common/units.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cmath>
#include <cstddef>

/**
 * @file case_config.hpp
 * @brief 五组分 PR、Søreide–Whitson 与 CPA 共同公开数据标定后的流动对比算例。
 *
 * 三次运行共用网格、岩石、井、初始 P-T-z 和输运离散，但每个 EOS 使用其
 * 纯组分与缔合参数仍保持各模型的公开参数化；可识别的六个二元参数对使用
 * 同一套 NIST ThermoML 训练数据分别标定，留出数据和流动结果均不参与拟合。
 * 未被公共数据覆盖的组分对保持原公开先验值。
 */
namespace CaseConfig
{

inline constexpr char name[] = "five_component_eos_tuned_compare";
inline constexpr double bar = 1.0e5;
inline constexpr double mD = 9.869232667160130e-16;
inline constexpr double secondsPerDay = 86400.0;
inline constexpr double gasConstant = MPMC::units::gasConstant;

// ============================================================================
// 1. Model: only the new fully-compositional O/G/W model is active.
// ============================================================================
struct Model
{
    static constexpr int numberOfComponents = 5;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase;
};

// ============================================================================
// 2. Three-dimensional structured reservoir.
// ============================================================================
struct Grid
{
    inline static constexpr const char *meshDirectory = "";

    // 20 x 20 x 5 = 2000 cells; each cell is 50 x 30 x 10 m.
    static constexpr int nx = 20;
    static constexpr int ny = 20;
    static constexpr int nz = 5;
    static constexpr double lx = 1000.0; // m
    static constexpr double ly = 600.0;  // m
    static constexpr double lz = 50.0;   // m
};

/**
 * @brief 简化确定性砂岩：单一直线高渗通道、单层隔层和一个 4x4 窗口。
 */
struct Rock
{
    static constexpr int channelCenter(int i)
    {
        // Exact zero-based counterpart of the MRST preview geometry.
        return static_cast<int>(4.5 + 11.0 * static_cast<double>(i - 1) / 17.0);
    }

    static constexpr bool inChannel(int i, int j, int k)
    {
        const int center = channelCenter(i);
        const int dj = j > center ? j - center : center - j;
        return k != 2 && dj <= 1;
    }

    static constexpr bool inBaffleWindow(int i, int j, int k)
    {
        return k == 2 && i >= 8 && i <= 11 && j >= 8 && j <= 11;
    }

    static constexpr double kx(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 150.0 * mD;
        if (k == 2)
            return 2.0 * mD;
        return (inChannel(i, j, k) ? 250.0 : 80.0) * mD;
    }
    static constexpr double ky(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 100.0 * mD;
        if (k == 2)
            return 1.0 * mD;
        return (inChannel(i, j, k) ? 160.0 : 50.0) * mD;
    }
    static constexpr double kz(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 10.0 * mD;
        if (k == 2)
            return 0.05 * mD;
        return (inChannel(i, j, k) ? 18.0 : 5.0) * mD;
    }
    static constexpr double porosity(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 0.22;
        if (k == 2)
            return 0.10;
        return inChannel(i, j, k) ? 0.24 : 0.19;
    }
};

// ============================================================================
// 3. Common component data and interaction baseline.
// ============================================================================
struct CommonFluid
{
    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int waterComponent = 0;
    static constexpr int co2Component = 1;
    static constexpr int methaneComponent = 2;
    static constexpr int ethaneComponent = 3;
    static constexpr int nButaneComponent = 4;

    inline static constexpr std::array<const char *, N> componentNames{
        "H2O", "CO2", "CH4", "C2H6", "nC4H10"};

    // Public pure-fluid reference data: CoolProp 8.0 fluid pages.
    // Units: K, Pa, m3/mol, dimensionless omega, kg/mol.
    inline static constexpr std::array<double, N> criticalTemperature{
        647.096, 304.128200003, 190.564002651, 305.322, 425.125};
    inline static constexpr std::array<double, N> criticalPressure{
        22.064000000e6, 7.377298373e6, 4.599200474e6, 4.872200000e6, 3.796000017e6};
    inline static constexpr std::array<double, N> criticalVolume{
        5.594803743e-5, 9.411848339e-5, 9.862771707e-5,
        1.458387816e-4, 2.549219298e-4};
    inline static constexpr std::array<double, N> acentricFactor{
        0.3442920843, 0.22394, 0.01142, 0.099, 0.200810094644};
    inline static constexpr std::array<double, N> molarMass{
        0.018015268, 0.0440098, 0.0160428, 0.03006904, 0.0581222};

    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr double eosU = 2.414213562373095;
    static constexpr double eosW = -0.414213562373095;

    static constexpr double temperature = 305.0; // K

    // Reference API values; fully-compositional reservoir properties are
    // evaluated from EOS state/composition at runtime.
    inline static constexpr std::array<double, 3> surfaceDensity{
        620.0, 1.80, 998.0};
    inline static constexpr std::array<double, 3> viscosity{
        3.0e-4, 1.5e-5, 8.0e-4};
    static constexpr double waterViscosity = 8.0e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;
};

struct PrFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;
    static constexpr int eosModelFlag = 5; // PR78 alpha branch

    // Isothermal volume translation anchored independently of the VLE fit.
    // At the case reference state (305 K, 6 MPa), c_H2O makes the pure-water
    // PR liquid volume equal the IAPWS-IF97 value rho=997.6771185 kg/m3.
    // Fugacity coefficients and phase-equilibrium conditions are unchanged.
    inline static constexpr std::array<double, N> componentVolumeTranslation{
        3.252997096253869e-6, 0.0, 0.0, 0.0, 0.0};

    // Common-data fit, 2026-08-25.  Unobserved H2O-C2/H2O-nC4, CO2-nC4 and
    // C2-nC4 pairs retain their previous public/reservoir priors.
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0,         -0.0391803910093841, -0.2000000000000000, 0.5000, 0.5000}},
        {{-0.0391803910093841, 0.0, 0.1100076406250000, 0.1743251250000000, 0.1277}},
        {{-0.2000000000000000, 0.1100076406250000, 0.0, 0.0022815210880000, 0.0108823211788800}},
        {{0.5000, 0.1743251250000000, 0.0022815210880000, 0.0, 0.0000}},
        {{0.5000, 0.1277, 0.0108823211788800, 0.0000, 0.0}}
    }};
};

/** @brief 含盐 SW：1 mol/kg NaCl，并采用 Chabab (2019) CO2 水相修正。 */
struct SwFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::SoreideWhitson;
    static constexpr int eosModelFlag = 1; // original Soreide-Whitson PR branch
    static constexpr double soreideWhitsonSalinityMolality = 1.0;

    // The SW water alpha gives a different unshifted liquid root from PR, so
    // its own H2O translation is anchored to the same independent IAPWS state.
    inline static constexpr std::array<double, N> componentVolumeTranslation{
        3.267331277345136e-6, 0.0, 0.0, 0.0, 0.0};

    // Water-role-only dissolved-CO2 correction.  The CO2 increment is obtained
    // at 305 K from Garcia's apparent molar volume
    // Vbar_CO2=35.32748593 cm3/mol, then expressed as an additional translation
    // relative to the already IAPWS-anchored SW aqueous volume at the common
    // initial aqueous composition (x_CO2=0.00302683).  It changes aqueous
    // density/volume but not SW alpha, BIP, Z or fugacity.
    inline static constexpr std::array<double, N>
        soreideWhitsonAqueousVolumeTranslation{
            0.0, 1.469715763209903e-6, 0.0, 0.0, 0.0};

    // Dry-pair values use the same common public training suite as PR/CPA.
    // Water-rich values remain SW correlations plus fitted constant offsets.
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0,    0.1896, 0.4850, 0.5000, 0.5000}},
        {{0.1896, 0.0,    0.1100076406250000, 0.1743251250000000, 0.1277}},
        {{0.4850, 0.1100076406250000, 0.0, 0.0022815210880000, 0.0108823211788800}},
        {{0.5000, 0.1743251250000000, 0.0022815210880000, 0.0, 0.0000}},
        {{0.5000, 0.1277, 0.0108823211788800, 0.0000, 0.0}}
    }};

    static constexpr double co2AqueousBipOffset = -0.0002375832200000;
    static constexpr double methaneAqueousBipOffset = -0.0378357132700000;

    static double soreideWhitsonAqueousWaterBip(
        int component,
        double temperature,
        double salinityMolality)
    {
        if (component == waterComponent)
            return 0.0;
        if (component == co2Component)
            return MPMC::SoreideWhitsonCorrelations::co2AqueousBipChabab2019(
                temperature, salinityMolality) + co2AqueousBipOffset;
        const double baseline = MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
            temperature,
            criticalTemperature[static_cast<std::size_t>(component)],
            acentricFactor[static_cast<std::size_t>(component)],
            salinityMolality);
        return baseline + (component == methaneComponent ? methaneAqueousBipOffset : 0.0);
    }
};

/**
 * @brief SRK-CPA：公开纯组分参数、B2 CO2 溶剂化和二元拟合参数。
 *
 * 来源：Qvistgaard et al., Fluid Phase Equilibria 570 (2023) 113796,
 * DOI 10.1016/j.fluid.2023.113796, Table 1，给出 b、Gamma、c1；本实现
 * 使用 a0=b*R*Gamma。二元参数和 B2（一负位点 CO2、实验交叉缔合能）来自
 * Tsivintzelis & Kontogeorgis, J. Supercrit. Fluids 104 (2015) 29-39,
 * DOI 10.1016/j.supflu.2015.05.015, Tables 2-3。该文在比较方案中推荐显式
 * CO2-water 溶剂化；这里逐项使用其 B2 参数，不对五组分结果再拟合。
 */
struct CpaFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::CubicPlusAssociation;
    static constexpr int eosModelFlag = 5;

    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0, 0.1584126867280001, 0.0204118168704000, 0.04415, 0.08750}},
        {{0.1584126867280001, 0.0, 0.0096233516000000, 0.0677109300736000, 0.11220}},
        {{0.0204118168704000, 0.0096233516000000, 0.0, -0.0045994250000000, -0.0218917961600000}},
        {{0.04415, 0.0677109300736000, -0.0045994250000000, 0.0, 0.00000}},
        {{0.08750, 0.11220, -0.0218917961600000, 0.00000, 0.0}}
    }};

    inline static constexpr std::array<double, N> cpaB{
        14.52e-6, 27.2e-6, 29.1e-6, 42.9e-6, 72.081e-6};
    inline static constexpr std::array<double, N> cpaGamma{
        1017.3, 1551.22, 959.02, 1544.54, 2193.08};
    inline static constexpr std::array<double, N> cpaA0 = [] {
        std::array<double, N> a{};
        for (std::size_t c = 0; c < N; ++c)
            a[c] = cpaB[c] * gasConstant * cpaGamma[c];
        return a;
    }();
    inline static constexpr std::array<double, N> cpaC1{
        0.6736, 0.7602, 0.4472, 0.5846, 0.7077};
    inline static constexpr std::array<double, N> cpaAssociationEnergy{
        2003.2 * gasConstant, 14200.0, 0.0, 0.0, 0.0};
    inline static constexpr std::array<double, N> cpaAssociationVolume{
        0.0692, 0.0162, 0.0, 0.0, 0.0};
    inline static constexpr std::array<int, N> cpaDonorSites{2, 0, 0, 0, 0};
    inline static constexpr std::array<int, N> cpaAcceptorSites{2, 1, 0, 0, 0};
    inline static constexpr std::array<std::array<double, N>, N>
        cpaCrossAssociationEnergy = [] {
            std::array<std::array<double, N>, N> values{};
            // 142.0 bar L/mol = 14 200 J/mol (B2 experimental value).
            values[waterComponent][co2Component] = 14200.0;
            return values;
        }();
    inline static constexpr std::array<std::array<double, N>, N>
        cpaCrossAssociationVolume = [] {
            std::array<std::array<double, N>, N> values{};
            values[waterComponent][co2Component] = 0.0162;
            return values;
        }();
    static constexpr auto cpaRadialDistribution =
        MPMC::CpaRadialDistribution::Simplified;
    static constexpr bool enableThermodynamicProfiler = true;
};

// Factory wrappers allow one executable to select PR/SW/CPA at runtime while
// retaining one common Runtime/Indices type.
struct PrFactoryConfig { using Fluid = PrFluid; };
struct SwFactoryConfig { using Fluid = SwFluid; };
struct CpaFactoryConfig { using Fluid = CpaFluid; };

// ============================================================================
// 4. Common P-T-z initialization. The preflight test requires all EOS to flash.
// ============================================================================
struct InitialState
{
    static constexpr double pressure = 60.0 * bar;
    static constexpr double temperature = CommonFluid::temperature;

    // [H2O, CO2, CH4, C2H6, nC4H10]. Wet condensable light-reservoir fluid.
    inline static constexpr std::array<double, CommonFluid::N> overallComposition{
        0.30, 0.10, 0.15, 0.15, 0.30};
};

// ============================================================================
// 5. Legacy optional physics are explicitly disabled by Model.
// ============================================================================
struct Dissolution
{
    static constexpr int component = CommonFluid::co2Component;
    static constexpr double waterMolarMass = 0.018015268;
    static constexpr double salinityMolality = 0.0;
    static constexpr double initialWaterCO2MoleFraction = 0.0;
};
struct Land { static constexpr double constant = 0.0; };
struct Adsorption
{
    static constexpr double rockDensity = 2650.0;
    inline static constexpr std::array<double, CommonFluid::N> thetaMax{};
    inline static constexpr std::array<double, CommonFluid::N> coefficient{};
    static constexpr double standardPressure = 101325.0;
    static constexpr double standardTemperature = 288.15;
};

struct Numerics
{
    static constexpr double fugacityScalingFactor = 1.0;
    static constexpr bool useVariableBounds = false;

    // 目标超算日志显示典型相边界平台：第一个 Newton 更新后残差明显下降，
    // 随后只在舍入误差量级变化，但旧流程仍把 SNES 跑满。保留较大的 max_it
    // 给真正困难但仍持续改善的求解，只对这种持续平台提前终止。
    static constexpr bool enableSnesStagnationGuard = true;
    static constexpr int snesStagnationMinimumIterations = 8;
    static constexpr int snesStagnationWindow = 5;
    static constexpr double snesStagnationRelativeImprovement = 1.0e-4;
};

// The 0.25-day grid defines common output targets for all EOS. Internal adaptive
// substeps are allowed between two targets; the transactional stepper rolls back
// a failed attempt and still lands exactly on every common output time.
struct Time
{
    static constexpr int numberOfSteps = 120;
    static constexpr double dtDays = 0.25;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = 1.0e-5;
    static constexpr double cutFactor = 0.5;
    static constexpr double growthFactor = 1.25;
    static constexpr double difficultShrinkFactor = 0.8;
    static constexpr int easyNonlinearIterations = 6;
    static constexpr int difficultNonlinearIterations = 14;
    static constexpr int maximumRetries = 14;
    static constexpr int maximumWellControlIterations = 6;
    static constexpr bool printAdaptiveSteps = true;
};

struct Output
{
    // If -result_dir is not supplied, the executable changes this to
    // ./results/pr, ./results/sw or ./results/cpa after parsing -eos.
    inline static constexpr const char *directory = "./results";
    static constexpr std::size_t every = 1;
    static constexpr bool printWells = true;
    static constexpr bool printWellPhaseDetails = true;
    static constexpr bool writeWellHistory = true;
    static constexpr bool printInventory = true;
    static constexpr bool enableComponentMassBalance = true;
    static constexpr bool printComponentMassBalance = true;
    static constexpr bool writeComponentMassBalance = true;
    static constexpr bool writeInventoryHistory = true;
    static constexpr bool writeMassTotals = true;
    static constexpr bool writeSolutionSnapshots = true;
    static constexpr bool writePhaseStateSnapshots = true;
    static constexpr bool printNewtonIterations = true;
    static constexpr bool writeSolverHistory = true;
    static constexpr bool printWellControlSwitches = true;
    static constexpr bool writeWellControlSwitchHistory = true;
    static constexpr bool printReservoirDiagnostics = true;
    static constexpr bool writeReservoirDiagnostics = true;
    static constexpr bool printFinalSummary = true;
    static constexpr bool writeFinalSummary = true;
};

struct Config final
{
    inline static constexpr const char *name = CaseConfig::name;
    using Model = CaseConfig::Model;
    using Grid = CaseConfig::Grid;
    using Fluid = CaseConfig::PrFluid; // common metadata/CSV naming baseline
    using InitialState = CaseConfig::InitialState;
    using Dissolution = CaseConfig::Dissolution;
    using Land = CaseConfig::Land;
    using Adsorption = CaseConfig::Adsorption;
    using Numerics = CaseConfig::Numerics;
    using Time = CaseConfig::Time;
    using Output = CaseConfig::Output;
    using Rock = CaseConfig::Rock;
};

} // namespace CaseConfig
