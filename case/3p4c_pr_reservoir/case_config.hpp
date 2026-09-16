#pragma once

#include <indices/model_config.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

/**
 * @file case_config.hpp
 * @brief H2O–CO2–CH4–nC16 三相四组分 PR 储层算例的网格、流体、时间步和求解参数。
 *
 * This is the first case in the project that uses the full three-phase model in
 * the actual PETSc reservoir path.  Initialization deliberately contains only
 * pressure, temperature and overall mole fractions.  Oil/gas/water-rich phase
 * fractions, phase compositions, Z factors and pore-volume saturations are all
 * produced by the PR P-T-z flash.
 *
 * The fluid set uses real component identities and a literature-style PR BIP
 * pattern for CO2/water/hydrocarbon systems.  It is a numerical/thermodynamic
 * benchmark, not a history-matched field-fluid characterization.
 */
namespace CaseConfig
{

inline constexpr char name[] = "3p4c_pr_reservoir";
inline constexpr double bar = 1.0e5;
inline constexpr double mD = 9.869232667160130e-16;
inline constexpr double secondsPerDay = 86400.0;

// ============================================================================
// 1. Model
// ============================================================================
struct Model
{
    static constexpr int numberOfComponents = 4;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase;
};

// ============================================================================
// 2. Structured reservoir grid and rock
// ============================================================================
struct Grid
{
    inline static constexpr const char *meshDirectory = "";

    // 30 x 15 x 1 = 450 cells.  The case is intentionally small enough for a
    // first PETSc smoke test while still having a 2-D displacement path.
    static constexpr int nx = 30;
    static constexpr int ny = 15;
    static constexpr int nz = 1;
    static constexpr double lx = 1200.0; // m
    static constexpr double ly = 600.0;  // m
    static constexpr double lz = 12.0;   // m
};

struct Rock
{
    static constexpr double kx = 100.0 * mD;
    static constexpr double ky = 100.0 * mD;
    static constexpr double kz = 10.0 * mD;
    static constexpr double porosity = 0.20;
};

// ============================================================================
// 3. Fully compositional PR fluid
// ============================================================================
struct Fluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;

    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int waterComponent = 0;
    static constexpr int co2Component = 1;
    static constexpr int methaneComponent = 2;
    static constexpr int heavyComponent = 3;

    inline static constexpr std::array<const char *, N> componentNames{
        "H2O", "CO2", "CH4", "nC16"};

    // Pure-component data used by the benchmark.  Units: K, Pa, m3/mol,
    // dimensionless omega, kg/mol.
    inline static constexpr std::array<double, N> criticalTemperature{
        647.30, 304.20, 190.564, 723.0};
    inline static constexpr std::array<double, N> criticalPressure{
        220.48e5, 73.76e5, 46.00e5, 14.10e5};
    inline static constexpr std::array<double, N> criticalVolume{
        5.6e-5, 9.4e-5, 9.9e-5, 9.0e-4};
    inline static constexpr std::array<double, N> acentricFactor{
        0.344, 0.225, 0.008, 0.742};
    inline static constexpr std::array<double, N> molarMass{
        0.01801528, 0.0440098, 0.016043, 0.2264412};

    // Symmetric PR binary-interaction matrix.  The water/CO2/hydrocarbon and
    // CO2/hydrocarbon magnitudes follow the interaction pattern used by the
    // Pang NWE full-three-phase benchmark; nC16 is used here as the heavy
    // hydrocarbon representative.
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0,    0.1896, 0.4850, 0.5000}},
        {{0.1896, 0.0,    0.1200, 0.0900}},
        {{0.4850, 0.1200, 0.0,    0.0000}},
        {{0.5000, 0.0900, 0.0000, 0.0}}
    }};

    // Peng-Robinson constants.  eosModelFlag=5 retains the PR78 kappa branch
    // for the high-acentric-factor heavy component while remaining the same PR
    // cubic EOS/mixing-rule backend used by the new full three-phase model.
    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr int eosModelFlag = 5;
    static constexpr double eosU = 2.414213562373095;
    static constexpr double eosW = -0.414213562373095;

    static constexpr double temperature = 350.0; // K, isothermal reservoir

    // Required FluidSystem reference values.  The fully-compositional path
    // evaluates reservoir density/viscosity from EOS composition and Z; these
    // remain useful surface/reference values and keep the common API complete.
    inline static constexpr std::array<double, 3> surfaceDensity{
        750.0, 100.0, 1000.0};
    inline static constexpr std::array<double, 3> viscosity{
        5.0e-4, 2.0e-5, 1.0e-4};
    static constexpr double waterViscosity = 1.0e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;
};

// ============================================================================
// 4. Initialization: P-T-z only
// ============================================================================
struct InitialState
{
    static constexpr double pressure = 200.0 * bar;
    static constexpr double temperature = Fluid::temperature;

    // Overall mole fractions [H2O, CO2, CH4, nC16].  The PR flash at 200 bar
    // and 350 K produces three active phases with approximately
    // So=0.37787, Sg=0.31743, Sw=0.30470.  Those saturations are diagnostics,
    // NOT user input.
    inline static constexpr std::array<double, Fluid::N> overallComposition{
        0.65, 0.30, 0.02, 0.03};
};

// ============================================================================
// 5. Disabled optional legacy physics (kept for common case_support API)
// ============================================================================
struct Dissolution
{
    static constexpr int component = Fluid::co2Component;
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

// ============================================================================
// 6. Short first reservoir test
// ============================================================================
struct Time
{
    static constexpr int numberOfSteps = 20;
    static constexpr double dtDays = 0.25;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = 1.0e-5;
    static constexpr double cutFactor = 0.5;
    static constexpr double growthFactor = 1.5;
    static constexpr double difficultShrinkFactor = 0.8;
    static constexpr int easyNonlinearIterations = 6;
    static constexpr int difficultNonlinearIterations = 14;
    static constexpr int maximumRetries = 12;
    static constexpr int maximumWellControlIterations = 6;
    static constexpr bool printAdaptiveSteps = true;
};

struct Output
{
    inline static constexpr const char *directory = "./results";
    static constexpr std::size_t every = 5;
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
