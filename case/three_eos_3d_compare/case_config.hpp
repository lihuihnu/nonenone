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
 * @brief PR、Søreide–Whitson 与 CPA 三 EOS 共用的三维规则网格对比算例参数。
 *
 * 算例只启用新的全组分 O/G/W 三相互溶模型。吸附、Land 气体滞留以及旧的
 * “独立水相 CO2 溶解变量”全部关闭。三个 EOS 共用完全相同的网格、岩石、井、
 * 初始 P-T-z、相对渗透率和非水相 BIP 基线；差异仅来自热力学后端。
 */
namespace CaseConfig
{

inline constexpr char name[] = "three_eos_3d_compare";
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

    // 36 x 24 x 8 = 6912 cells; cell sizes are about 33.3 x 30 x 12 m.
    static constexpr int nx = 36;
    static constexpr int ny = 24;
    static constexpr int nz = 8;
    static constexpr double lx = 1200.0; // m
    static constexpr double ly = 720.0;  // m
    static constexpr double lz = 96.0;   // m
};

/**
 * @brief 确定性分层砂岩模型，含一层低渗隔夹层和一条斜向高渗通道。
 *
 * 物理：主层 60--220 mD，垂向渗透率约为水平渗透率的 5--10%；第 4 层
 * 为低渗隔夹层，但在模型中部保留有限导流窗口。斜向高渗通道连接注采井附近区域，
 * 下部注入流体需先向窗口汇聚再进入上部生产层，从而形成可重复的三维流线聚焦。
 */
struct Rock
{
    inline static constexpr std::array<double, Grid::nz> layerKxMd{
        220.0, 190.0, 150.0, 5.0, 95.0, 140.0, 110.0, 70.0};
    inline static constexpr std::array<double, Grid::nz> layerKyMd{
        180.0, 160.0, 125.0, 4.0, 80.0, 115.0, 90.0, 60.0};
    inline static constexpr std::array<double, Grid::nz> layerKzMd{
        18.0, 15.0, 12.0, 0.45, 7.5, 11.0, 8.5, 5.0};
    inline static constexpr std::array<double, Grid::nz> layerPorosity{
        0.225, 0.215, 0.205, 0.120, 0.180, 0.200, 0.190, 0.175};

    static constexpr int channelCenter(int i)
    {
        return 4 + (15 * i) / (Grid::nx - 1);
    }

    static constexpr bool inChannel(int i, int j, int k)
    {
        const int center = channelCenter(i);
        const int dj = j > center ? j - center : center - j;
        return k != 3 && dj <= 2;
    }

    static constexpr bool inBaffleWindow(int i, int j, int k)
    {
        if (k != 3 || i < 16 || i > 20)
            return false;
        const int center = channelCenter(i);
        const int dj = j > center ? j - center : center - j;
        return dj <= 2;
    }

    static constexpr double kx(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 80.0 * mD;
        return layerKxMd[static_cast<std::size_t>(k)] *
               (inChannel(i, j, k) ? 2.0 : 1.0) * mD;
    }
    static constexpr double ky(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 65.0 * mD;
        return layerKyMd[static_cast<std::size_t>(k)] *
               (inChannel(i, j, k) ? 2.0 : 1.0) * mD;
    }
    static constexpr double kz(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 6.0 * mD;
        return layerKzMd[static_cast<std::size_t>(k)] *
               (inChannel(i, j, k) ? 1.5 : 1.0) * mD;
    }
    static constexpr double porosity(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 0.18;
        const double base = layerPorosity[static_cast<std::size_t>(k)];
        return base + (inChannel(i, j, k) ? 0.015 : 0.0);
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

    // Common non-aqueous BIP baseline used by all three EOS so the comparison
    // does not silently retune each backend. Water interactions follow the
    // established SW non-aqueous values; the light-hydrocarbon/CO2 block uses
    // the SPE3-style pattern already validated in this project. Unlisted pairs
    // are deliberately zero rather than fitted to this synthetic reservoir.
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0,    0.1896, 0.4850, 0.5000, 0.5000}},
        {{0.1896, 0.0,    0.1000, 0.1300, 0.1277}},
        {{0.4850, 0.1000, 0.0,    0.0000, 0.09281}},
        {{0.5000, 0.1300, 0.0000, 0.0,    0.0000}},
        {{0.5000, 0.1277, 0.09281,0.0000, 0.0}}
    }};

    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr int eosModelFlag = 1; // original PR76 kappa branch
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
};

/** @brief Søreide–Whitson PR：采用淡水 alpha 与水相 H2O–溶质 BIP 关联式。 */
struct SwFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::SoreideWhitson;
    static constexpr double soreideWhitsonSalinityMolality = 0.0;

    static double soreideWhitsonAqueousWaterBip(
        int component,
        double temperature,
        double salinityMolality)
    {
        if (component == waterComponent)
            return 0.0;
        if (component == co2Component)
            return MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
                temperature,
                criticalTemperature[static_cast<std::size_t>(component)],
                salinityMolality);
        return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
            temperature,
            criticalTemperature[static_cast<std::size_t>(component)],
            acentricFactor[static_cast<std::size_t>(component)],
            salinityMolality);
    }
};

/**
 * @brief SRK-CPA：显式采用公开 CPA 纯组分表和 4C-water 缔合参数。
 *
 * 来源：Qvistgaard et al., Fluid Phase Equilibria 570 (2023) 113796,
 * DOI 10.1016/j.fluid.2023.113796, Table 1。表中给 b、Gamma、c1；
 * 本实现按该 CPA 参数定义使用 a0=b*R*Gamma。只有水具有 4C 缔合位点。
 */
struct CpaFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::CubicPlusAssociation;

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
        2003.2 * gasConstant, 0.0, 0.0, 0.0, 0.0};
    inline static constexpr std::array<double, N> cpaAssociationVolume{
        0.0692, 0.0, 0.0, 0.0, 0.0};
    inline static constexpr std::array<int, N> cpaDonorSites{2, 0, 0, 0, 0};
    inline static constexpr std::array<int, N> cpaAcceptorSites{2, 0, 0, 0, 0};
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
// 4. Common P-T-z initialization. All three EOS resolve O/G/W at this state.
// ============================================================================
struct InitialState
{
    static constexpr double pressure = 60.0 * bar;
    static constexpr double temperature = CommonFluid::temperature;

    // [H2O, CO2, CH4, C2H6, nC4H10]. This is an intentionally wet,
    // condensable light-reservoir mixture used to exercise all O/G/W phases.
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

// The 0.5-day grid defines common output targets for all EOS. Internal adaptive
// substeps are allowed between two targets; the transactional stepper rolls back
// a failed attempt and still lands exactly on every common output time.
struct Time
{
    static constexpr int numberOfSteps = 120;
    static constexpr double dtDays = 0.5;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = 1.0e-5;
    static constexpr double cutFactor = 0.5;
    static constexpr double growthFactor = 1.25;
    static constexpr double difficultShrinkFactor = 0.8;
    static constexpr int easyNonlinearIterations = 6;
    static constexpr int difficultNonlinearIterations = 14;
    static constexpr int maximumRetries = 12;
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
