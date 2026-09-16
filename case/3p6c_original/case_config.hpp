#pragma once

#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

/**
 * @file case_config.hpp
 * @brief 三相六组分原始对照算例的网格、流体、时间步和求解参数。
 *
 * 新算例最常改本文件：Model / Grid / Rock / Fluid / InitialState / Time / Output。
 * 井参数只放 well_config.hpp，main 只保留调用流程。
 */
namespace CaseConfig
{

inline constexpr char name[] = "3p6c_original";
inline constexpr double bar = 1.0e5;
inline constexpr double mD = 9.869232667160130e-16;

// ============================================================================
// 1. 模型开关：新增物理功能时优先从这里开关
// ============================================================================
struct Model
{
    static constexpr int numberOfComponents = 6;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
};

// ============================================================================
// 2. 结构网格与岩石
// ============================================================================
struct Grid
{
    // StructuredGrid 不读取外部网格目录；保留统一字段让 case_support API 更直白。
    inline static constexpr const char *meshDirectory = "";

    static constexpr int nx = 20;
    static constexpr int ny = 20;
    static constexpr int nz = 2;
    static constexpr double lx = 1000.0; // m
    static constexpr double ly = 1000.0; // m
    static constexpr double lz = 1.0;    // m
};

struct Rock
{
    static constexpr double kx = 50.0 * mD;
    static constexpr double ky = 50.0 * mD;
    static constexpr double kz = 50.0 * mD;
    static constexpr double porosity = 0.25;
};

// ============================================================================
// 3. 组分、EOS、流体
// ============================================================================
struct Fluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;

    static constexpr std::size_t N = Model::numberOfComponents;

    inline static constexpr std::array<const char *, N> componentNames{
        "N2/CH4", "CO2", "C2-5", "C6-13", "C14-24", "C25-80"};

    inline static constexpr std::array<double, N> criticalTemperature{
        189.515, 304.2, 387.607, 597.497, 698.515, 875.0};
    inline static constexpr std::array<double, N> criticalPressure{
        4580011.59, 7386592.50, 4095515.97,
        3345244.875, 1768374.5625, 1169006.79};
    inline static constexpr std::array<double, N> criticalVolume{
        9.97012032965401e-5, 9.26344713533338e-5,
        2.17076707259486e-4, 3.81162235869935e-4,
        7.21410148317871e-4, 1.13570073874421e-3};
    inline static constexpr std::array<double, N> acentricFactor{
        0.00854, 0.228, 0.16733, 0.38609, 0.80784, 1.23141};
    inline static constexpr std::array<double, N> molarMass{
        0.0161594, 0.04401, 0.0455725, 0.11774, 0.248827, 0.48152};

    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0,        0.00070981, 0.00077754, 0.0100, 0.0110, 0.0110}},
        {{0.00070981, 0.0,        0.1500,     0.1500, 0.1500, 0.1500}},
        {{0.00077754, 0.1500,     0.0,        0.0,    0.0,    0.0}},
        {{0.0100,     0.1500,     0.0,        0.0,    0.0,    0.0}},
        {{0.0110,     0.1500,     0.0,        0.0,    0.0,    0.0}},
        {{0.0110,     0.1500,     0.0,        0.0,    0.0,    0.0}}
    }};

    static constexpr double eosOmegaA = 0.4572355;
    static constexpr double eosOmegaB = 0.0779691;
    static constexpr int eosModelFlag = 1;
    static constexpr double eosU = 2.414213562373095;
    static constexpr double eosW = -0.414213562373095;

    static constexpr double temperature = 387.45; // K
    inline static constexpr std::array<double, 3> surfaceDensity{800.0, 1.86906, 1000.0};
    inline static constexpr std::array<double, 3> viscosity{0.5, 0.7, 2.0e-4};
    static constexpr double waterViscosity = 2.0e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;
};

// ============================================================================
// 4. 初始状态：严格沿用旧 3p6c_verify.hpp
// ============================================================================
struct InitialState
{
    static constexpr double pressure = 150.0 * bar;
    // 以下数值逐项来自原 MPMC_LH/3p6c_verify.hpp。
    static constexpr double waterSaturation = 0.2;
    static constexpr double oilSaturation = 0.5308359;
    static constexpr double gasSaturation = 0.2691692;

    inline static constexpr std::array<double, Fluid::N> oilComposition{
        0.3246914, 0.0128351, 0.2278401,
        0.2606985, 0.1134144, 0.0605206};
    inline static constexpr std::array<double, Fluid::N> gasComposition{
        0.808671, 0.025284, 0.1487798,
        0.0175878, 0.0006759, 1.51e-06};

    // 原文件中的初值存在轻微非闭合；为做数值回归，按原值保留，不在 case 层归一化。
    static constexpr bool preserveReferenceValues = true;

    // true：保留原算例直接写入 phaseState 的 K/z/L；false：重新 flash。
    static constexpr bool useLegacySecondaryState = true;
    static constexpr int hydrocarbonPhaseFlag = 0; // 0 = oil+gas two-phase
    inline static constexpr std::array<double, Fluid::N> equilibriumRatio{
        0.401512, 0.507637, 1.531391,
        14.82269, 167.797603, 4.00798e-8};
    inline static constexpr std::array<double, Fluid::N> overallComposition{
        0.463, 0.0164, 0.2052, 0.19108, 0.08113, 0.04319};
    static constexpr double liquidMoleFraction = 0.727863364726590;
};

// ============================================================================
// 5. 可选物理参数：即使功能关闭也保留，未来只需打开 Model 开关
// ============================================================================
struct Dissolution
{
    static constexpr int component = 1; // CO2
    static constexpr double waterMolarMass = 0.01801528;
    static constexpr double salinityMolality = 0.0;
    static constexpr double initialWaterCO2MoleFraction = 0.0;
};

struct Land
{
    static constexpr double constant = 0.4;
};

struct Adsorption
{
    static constexpr double rockDensity = 2650.0; // kg/m3
    inline static constexpr std::array<double, Fluid::N> thetaMax{};
    inline static constexpr std::array<double, Fluid::N> coefficient{};
    static constexpr double standardPressure = 101325.0;
    static constexpr double standardTemperature = 288.15;
};

struct Numerics
{
    static constexpr double fugacityScalingFactor = 1.0;
    static constexpr bool useVariableBounds = false;
};

// ============================================================================
// 6. 时间和输出：可由 run.sh 环境变量/命令行临时覆盖
// ============================================================================
struct Time
{
    static constexpr int numberOfSteps = 100;
    // 原 run.sh: -numSteps 100 -dt 10。AdaptiveTimeStepper 只在需要时切内部子步。
    static constexpr double dtDays = 10.0;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = dtDays / 1024.0;
    static constexpr double cutFactor = 0.5;
    static constexpr double growthFactor = 2.0;
    static constexpr double difficultShrinkFactor = 0.8;
    static constexpr int easyNonlinearIterations = 10;
    static constexpr int difficultNonlinearIterations = 20;
    static constexpr int maximumRetries = 10;
    static constexpr int maximumWellControlIterations = 6;
    static constexpr bool printAdaptiveSteps = true;
};

struct Output
{
    inline static constexpr const char *directory = "./results";
    static constexpr std::size_t every = 25;
    static constexpr bool printWells = true;
    static constexpr bool printWellPhaseDetails = true;
    static constexpr bool writeWellHistory = true;
    static constexpr bool printInventory = true;
    static constexpr bool enableComponentMassBalance = false;
    static constexpr bool printComponentMassBalance = false;
    static constexpr bool writeComponentMassBalance = false;
    static constexpr bool writeInventoryHistory = true;
    static constexpr bool writeMassTotals = true;
    static constexpr bool writeSolutionSnapshots = true;
    static constexpr bool writePhaseStateSnapshots = true;

    // Solver diagnostics: custom MPMC monitor replaces noisy PETSc default monitors.
    static constexpr bool printNewtonIterations = true;
    static constexpr bool writeSolverHistory = true;

    // Rare but important well-control transitions get their own screen/CSV log.
    static constexpr bool printWellControlSwitches = true;
    static constexpr bool writeWellControlSwitchHistory = true;

    // Reservoir-wide state diagnostics are printed at OUT_STEP and written every output time.
    static constexpr bool printReservoirDiagnostics = true;
    static constexpr bool writeReservoirDiagnostics = true;

    // End-of-run statistics for paper tables and performance analysis.
    static constexpr bool printFinalSummary = true;
    static constexpr bool writeFinalSummary = true;
};


/** @brief 供 case_support.hpp 使用的类型化配置视图。 */
struct Config final
{
    inline static constexpr const char *name = CaseConfig::name;
    using Model = CaseConfig::Model;
    using Grid = CaseConfig::Grid;
    using Fluid = CaseConfig::Fluid;
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
