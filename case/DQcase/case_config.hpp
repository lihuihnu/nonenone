#pragma once

#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

/**
 * @file case_config.hpp
 * @brief DQ 网格三相多组分超算算例的网格、流体、时间步和求解参数。
 *
 * 网格/岩石来自 DQ_data；井只在 well_config.hpp 中配置。
 * InitialState 逐项对齐用户提供的原 MPMC_LH/DQcase.hpp。
 */
namespace CaseConfig
{

inline constexpr char name[] = "DQcase";
inline constexpr double bar = 1.0e5;
inline constexpr double secondsPerDay = 86400.0;

struct Model
{
    static constexpr int numberOfComponents = 8;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = true;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = true;
};

struct Grid
{
    // Machine-specific DQ_data paths live in config/hpc.mk or hpc.local.mk.
    // The managed run.sh passes that profile value through -mesh_dir. Direct
    // executable runs may also provide -mesh_dir or -grdecl explicitly.
    inline static constexpr const char *meshDirectory = "";
    static constexpr bool replaceZeroPermeabilityWithMean = false;
    static constexpr bool replaceZeroPorosityWithMean = false;
};

struct Fluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;

    static constexpr std::size_t N = Model::numberOfComponents;
    inline static constexpr std::array<const char *, N> componentNames{
        "CO2", "N2+C1", "C2-nC4", "iC5+nC5+",
        "C7-C17", "C18-C22", "C23-C27", "C28-C80"};

    inline static constexpr std::array<double, N> criticalTemperature{
        304.2000, 183.8634, 372.9258, 488.3443,
        654.7271, 736.7347, 773.6587, 828.3198};
    inline static constexpr std::array<double, N> criticalPressure{
        7376460.000000, 4474024.626750, 4194267.315000, 3153252.238500,
        1994141.861250, 1481527.540500, 1357534.111500, 1237614.960750};
    inline static constexpr std::array<double, N> criticalVolume{
        9.257314158824152e-05, 9.908432031274359e-05,
        2.069828561549373e-04, 3.476491681729444e-04,
        6.824240059315388e-04, 9.509096452966007e-04,
        9.950112186761062e-04, 1.001601617958e-03};
    inline static constexpr std::array<double, N> acentricFactor{
        0.2250000, 0.0113474, 0.1522270, 0.2703830,
        0.5831180, 0.8456000, 0.9733940, 1.0845000};
    inline static constexpr std::array<double, N> molarMass{
        0.04400980, 0.01679355, 0.04357611, 0.07925293,
        0.16437050, 0.27655590, 0.34112370, 0.44980100};

    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0,       0.1114088, 0.1200000, 0.1200000, 0.1000000, 0.1000000, 0.1000000, 0.09999999}},
        {{0.1114088, 0.0,       0.004812325, 0.005515154, 0.005016739, 0.005016738, 0.005016738, 0.005016736}},
        {{0.1200000, 0.004812325, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{0.1200000, 0.005515154, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{0.1000000, 0.005016739, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{0.1000000, 0.005016738, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{0.1000000, 0.005016738, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{0.09999999,0.005016736, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}
    }};

    static constexpr double eosOmegaA = 0.4572355;
    static constexpr double eosOmegaB = 0.0779691;
    static constexpr int eosModelFlag = 1;
    static constexpr double eosU = 2.414213562373095;
    static constexpr double eosW = -0.414213562373095;

    static constexpr double temperature = 371.7611; // K

    // 当前 DQ 对齐基线；相顺序 Oil/Gas/Water。
    inline static constexpr std::array<double, 3> surfaceDensity{857.16, 1.86894, 1000.0};
    inline static constexpr std::array<double, 3> viscosity{5.0e-4, 7.0e-4, 1.0e-4};
    static constexpr double waterViscosity = 1.0e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;
};

struct InitialState
{
    // 原 DQcase.hpp 的 initDQcase() 初值。
    static constexpr double pressure = 22914370.0; // Pa
    static constexpr double waterSaturation = 0.55;
    static constexpr double oilSaturation = 0.45;
    static constexpr double gasSaturation = 0.0;

    inline static constexpr std::array<double, Fluid::N> oilComposition{
        1.0e-6, 0.2131767, 0.0440212, 0.0434585,
        0.3071616, 0.1744741, 0.1251466, 0.0925433};
    inline static constexpr std::array<double, Fluid::N> gasComposition{
        1.6201622674232575e-06, 9.8684419813110624e-01,
        1.2020673609374258e-02, 1.0696525651381764e-03,
        6.2713864628792385e-05, 1.0149916395226942e-06,
        1.1859652422041062e-07, 8.0793214561843329e-09};

    // 原文件小数截断后并非严格闭合；做基准回归时保留原始输入。
    static constexpr bool preserveReferenceValues = true;
    static constexpr bool useLegacySecondaryState = true;
    static constexpr int hydrocarbonPhaseFlag = 1; // 1 = liquid-only
    static constexpr double liquidMoleFraction = 1.0;
    inline static constexpr std::array<double, Fluid::N> equilibriumRatio{
        1.0646752108436826e+00, 3.0434242955373900e+00,
        1.7952313053890312e-01, 1.6181595181032193e-02,
        1.3422756169463679e-04, 3.8245926930631031e-06,
        6.2302621228421948e-07, 5.7395641458852598e-08};
    inline static constexpr std::array<double, Fluid::N> overallComposition{
        1.0004589895361992e-06, 2.1317879777036294e-01,
        4.4021649539577566e-02, 4.3458959545462736e-02,
        3.0717069678730169e-01, 1.7447579817515760e-01,
        1.2514789869107812e-01, 9.2545199032069744e-02};
};

struct Dissolution
{
    static constexpr int component = 0; // CO2
    static constexpr double waterMolarMass = 0.01801528;
    static constexpr double salinityMolality = 0.0;
    static constexpr double initialWaterCO2MoleFraction = 3.373e-8;
};

struct Land
{
    static constexpr double constant = 0.4;
};

struct Adsorption
{
    static constexpr double rockDensity = 2650.0;
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

struct Time
{
    static constexpr int numberOfSteps = 365;
    // 原 DQcase_adaptive run_sw_0.5.sh: -numSteps 365 -dt 10。
    static constexpr double dtDays = 10.0;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = 1.0e-9;
    static constexpr double cutFactor = 0.5;
    static constexpr double growthFactor = 1.5;
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
    static constexpr std::size_t every = 60;
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
};

} // namespace CaseConfig
