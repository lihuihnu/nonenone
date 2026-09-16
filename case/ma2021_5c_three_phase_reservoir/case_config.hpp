/**
 * @file case_config.hpp
 * @brief Ma 2021 五组分油气水三相文献算例的网格、流体、时间步和求解参数。
 */
#pragma once

#include <indices/model_config.hpp>
#include <common/units.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

/**
 * @file case_config.hpp
 * @brief Ma 等（2021）H2O-C1-C6-C10-C15 真实 O/G/W 互溶基准算例。
 *
 * Thermodynamic source:
 *   Ma et al., "Three-Phase Equilibrium Calculations of
 *   Water/Hydrocarbon/Nonhydrocarbon Systems Based on the Equation of State
 *   (EOS) in Thermal Processes", ACS Omega 6(50) (2021) 34406-34415,
 *   DOI: 10.1021/acsomega.1c04522, Tables 7-9.
 *
 * The paper's three-phase split is at 13.79 bar and 366.5 K for the overall
 * mixture H2O/C1/C6/C10/C15 = 0.10/0.10/0.20/0.40/0.20.  Unlike a
 * free-water/Henry shortcut, the paper places all components, including H2O,
 * in the common EOS phase-equilibrium problem and reports finite H2O in both
 * hydrocarbon phases and finite C1 in the aqueous phase.
 *
 * Important model-compatibility boundary:
 *   Ma et al. use a modified Peng-Robinson/Soreide-Whitson treatment with
 *   phase-dependent aqueous/non-aqueous BIPs and a modified H2O alpha term.
 *   MPMC_SCW v23 currently has one ordinary-PR BIP matrix shared by all phases.
 *   Therefore the P-T-z/component data and the NON-AQUEOUS BIPs below are
 *   literature values, while the resulting MPMC phase fractions are an
 *   ordinary-PR comparison, not a claim of exact Table-9 reproduction.
 *
 * The paper is a thermodynamic flash benchmark rather than a reservoir-flow
 * geometry benchmark.  Grid, rock and wells below are deliberately small MPMC
 * transport wrappers and are labelled as such.
 */
namespace CaseConfig
{

inline constexpr char name[] = "ma2021_5c_three_phase_reservoir";
inline constexpr double bar = 1.0e5;
inline constexpr double mD = 9.869232667160130e-16;
inline constexpr double secondsPerDay = 86400.0;
inline constexpr double gasConstant = MPMC::units::gasConstant;

// ============================================================================
// 1. Model: every component, including H2O, belongs to the common O/G/W EOS.
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
// 2. MPMC reservoir wrapper (NOT a parameter table from Ma et al.).
// ============================================================================
struct Grid
{
    inline static constexpr const char *meshDirectory = "";
    static constexpr int nx = 40;
    static constexpr int ny = 1;
    static constexpr int nz = 1;
    static constexpr double lx = 200.0; // m
    static constexpr double ly = 10.0;  // m
    static constexpr double lz = 10.0;  // m
};

struct Rock
{
    static constexpr double kx = 100.0 * mD;
    static constexpr double ky = 100.0 * mD;
    static constexpr double kz = 10.0 * mD;
    static constexpr double porosity = 0.20;
};

// ============================================================================
// 3. Ma et al. (2021), Tables 7-8: five-component three-phase fluid.
// ============================================================================
struct Fluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;

    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int waterComponent = 0;
    static constexpr int methaneComponent = 1;
    static constexpr int hexaneComponent = 2;
    static constexpr int decaneComponent = 3;
    static constexpr int pentadecaneComponent = 4;

    inline static constexpr std::array<const char *, N> componentNames{
        "H2O", "C1", "C6", "C10", "C15"};

    // Table 7. Units converted to SI where needed.
    inline static constexpr std::array<double, N> criticalTemperature{
        647.3, 190.6, 507.5, 622.1, 718.6};
    // The machine-readable PMC rendering prints C15 Pc as "08.49" in
    // Table 7.  This is a transcription defect: the same C15 property is
    // 18.49 bar in Table 4 immediately above, and 18.49 bar is used here.
    inline static constexpr std::array<double, N> criticalPressure{
        220.47e5, 46.00e5, 32.89e5, 25.34e5, 18.49e5};
    inline static constexpr std::array<double, N> acentricFactor{
        0.344, 0.008, 0.275, 0.444, 0.651};
    inline static constexpr std::array<double, N> molarMass{
        0.018, 0.016, 0.086, 0.134, 0.206};

    // Table 7 does not report critical volumes.  Natural's current LBC-style
    // viscosity closure requires Vc, although PR fugacity does not.  Use one
    // transparent project-side estimate Vc = Zc*R*Tc/Pc with Zc=0.27 rather
    // than silently importing values from another source.
    inline static constexpr std::array<double, N> criticalVolume{
        0.27 * gasConstant * criticalTemperature[0] / criticalPressure[0],
        0.27 * gasConstant * criticalTemperature[1] / criticalPressure[1],
        0.27 * gasConstant * criticalTemperature[2] / criticalPressure[2],
        0.27 * gasConstant * criticalTemperature[3] / criticalPressure[3],
        0.27 * gasConstant * criticalTemperature[4] / criticalPressure[4]};

    // Table 8: NON-AQUEOUS BIPs.  Ma et al. additionally use a second,
    // temperature-dependent aqueous BIP set; the current MPMC cubic-EOS API
    // has only one matrix, so this table is the explicit comparison choice.
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0,      0.4850,   0.4800,   0.4800,   0.4800}},
        {{0.4850,   0.0,      0.0,      0.0,      0.0}},
        {{0.4800,   0.0,      0.0,      0.002866, 0.010970}},
        {{0.4800,   0.0,      0.002866, 0.0,      0.002657}},
        {{0.4800,   0.0,      0.010970, 0.002657, 0.0}}
    }};

    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    // Use the ordinary PR76 alpha branch so this case remains an explicit
    // baseline against Ma et al.'s modified-PR treatment.
    static constexpr int eosModelFlag = 1;
    static constexpr double eosU = 2.414213562373095;
    static constexpr double eosW = -0.414213562373095;

    static constexpr double temperature = 366.5; // K, Table 9 flash condition

    // Surface/reference values required by the common well/output API. They are
    // project-side engineering values, not Table-7 data and do not enter PR
    // fugacity equilibrium.
    inline static constexpr std::array<double, 3> surfaceDensity{
        750.0, 0.72, 1000.0};
    inline static constexpr std::array<double, 3> viscosity{
        5.0e-4, 1.5e-5, 3.0e-4};
    static constexpr double waterViscosity = 3.0e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;
};

// ============================================================================
// 4. Literature P-T-z initial state, Ma et al. Tables 7-9.
// ============================================================================
struct InitialState
{
    static constexpr double pressure = 13.79 * bar;
    static constexpr double temperature = Fluid::temperature;
    inline static constexpr std::array<double, Fluid::N> overallComposition{
        0.10, 0.10, 0.20, 0.40, 0.20};
};

// Reference values from Ma et al. Table 9.  These are kept in the production
// case config so tests can compare the exact case input against the paper.
struct LiteratureReference
{
    // Order: Oil-rich, gas-rich, water-rich; components H2O,C1,C6,C10,C15.
    inline static constexpr std::array<std::array<double, Fluid::N>, 3> composition{{
        {{0.00622,  0.03961,  0.13372, 0.27329,  0.54717}},
        {{0.059951, 0.908347, 0.028802, 0.002770, 0.000119}},
        {{0.999587, 0.000004, 0.0,      0.0,      0.0}}
    }};
    inline static constexpr std::array<double, 3> phaseMoleFraction{
        0.73102, 0.07813, 0.19085};
};

// ============================================================================
// 5. Disabled optional legacy physics (required by common case API).
// ============================================================================
struct Dissolution
{
    static constexpr int component = Fluid::methaneComponent;
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
// 6. Small isothermal flow wrapper for PETSc/transport/mass-balance validation.
// ============================================================================
struct Time
{
    static constexpr int numberOfSteps = 20;
    static constexpr double dtDays = 0.10;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = 1.0e-6;
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
