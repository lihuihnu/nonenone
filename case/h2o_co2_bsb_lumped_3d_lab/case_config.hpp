/**
 * @file case_config.hpp
 * @brief 定义实验室尺度六组分超临界水三维流动算例的物理与数值参数。
 */
#pragma once

#define MPMC_H2O_CO2_NC10_INITIAL_PRESSURE_PA 28.0e6
#define MPMC_H2O_CO2_NC10_TEMPERATURE_K 653.15

#include "../h2o_co2_nc10_2d_benchmark/benchmark_common.hpp"

#include <indices/model_config.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

namespace LabScw3D
{

inline constexpr char name[] = "h2o_co2_bsb_lumped_3d_lab";
inline constexpr double secondsPerDay = 86400.0;

struct Model
{
    static constexpr int numberOfComponents = 6;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase;
};

// Reuse the H2O-CO2-nC10 comparison geometry and rock verbatim.  Keeping
// aliases here (instead of copying numbers) makes accidental grid drift a
// compile-time-visible source change in the common benchmark.
using Grid = BenchmarkCommon::Grid;
using Rock = BenchmarkCommon::Rock;

struct CommonFluid
{
    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int waterComponent = 0;
    static constexpr int co2Component = 1;

    inline static constexpr std::array<const char *, N> componentNames{
        "H2O", "CO2", "L_C1_C6", "M_C7_C15", "H_C16_C27", "XH_C28plus"};

    // H2O, CO2 and BSB pseudo-component data are taken from the PPT.  L is a
    // PPT-mole-weighted C1/C2-3/C4-6 lump.  Vc for BSB lumps uses the
    // preregistered Zc=0.27 estimate and must be calibrated before predictive use.
    inline static constexpr std::array<double, N> criticalTemperature{
        647.30, 304.20, 354.1916431226766, 605.78, 751.00, 942.50};
    inline static constexpr std::array<double, N> criticalPressure{
        22.048e6, 7.376e6, 4.065799256505576e6, 2.175e6, 1.654e6, 1.642e6};
    inline static constexpr std::array<double, N> criticalVolume{
        5.594803743e-5, 9.411848339e-5, 1.95564637471291e-4,
        6.252498825299838e-4, 1.0193008374141067e-3,
        1.2885644791440596e-3};
    inline static constexpr std::array<double, N> acentricFactor{
        0.344, 0.225, 0.1498936802973978, 0.618, 0.957, 1.268};
    inline static constexpr std::array<double, N> molarMass{
        0.018015, 0.044010, 0.04606110037174722, 0.140960, 0.280990, 0.519620};

    // Preregistered engineering baseline, not fitted to a flow result.
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0,    0.1896, 0.5000, 0.5000, 0.5000, 0.5000}},
        {{0.1896, 0.0,    0.12264602230483274, 0.0900, 0.0900, 0.0900}},
        {{0.5000, 0.12264602230483274, 0.0, 0.0, 0.0, 0.0}},
        {{0.5000, 0.0900, 0.0, 0.0, 0.0, 0.0}},
        {{0.5000, 0.0900, 0.0, 0.0, 0.0, 0.0}},
        {{0.5000, 0.0900, 0.0, 0.0, 0.0, 0.0}}
    }};

    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr int eosModelFlag = 1;
    static constexpr double eosU = 2.414213562373095;
    static constexpr double eosW = -0.414213562373095;
    static constexpr double temperature = BenchmarkCommon::temperature;

    // Reference densities only define the surface-rate representation.  The
    // experiment controls ReservoirTotalRate and evaluates reservoir properties
    // from the selected EOS at runtime.
    inline static constexpr std::array<double, 3> surfaceDensity{
        750.0, 1.8, 985.4040020947351};
    inline static constexpr std::array<double, 3> viscosity{
        3.0e-4, 1.86e-5, 4.69091e-4};
    static constexpr double waterViscosity = 4.69091e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;

    static constexpr bool useIapwsGarciaAqueousVolume = true;
    static constexpr double aqueousVolumeWaterMolarMass = 0.018015268;
    static constexpr double aqueousVolumeMaximumUnsupportedMoleFraction = 1.0e-4;
    static constexpr bool useMcBrideWrightAqueousViscosity = true;
    static constexpr double aqueousViscosityMaximumUnsupportedMoleFraction = 1.0e-4;
};

struct PrFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;
};

struct SwFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::SoreideWhitson;
    static constexpr double soreideWhitsonSalinityMolality = 0.0;

    static double soreideWhitsonAqueousWaterBip(
        int component, double temperature, double salinityMolality)
    {
        if (component == waterComponent)
            return 0.0;
        if (component == co2Component)
            return MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
                temperature, criticalTemperature[static_cast<std::size_t>(component)],
                salinityMolality);
        return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
            temperature, criticalTemperature[static_cast<std::size_t>(component)],
            acentricFactor[static_cast<std::size_t>(component)], salinityMolality);
    }
};

struct PrFactoryConfig { using Fluid = PrFluid; };
struct SwFactoryConfig { using Fluid = SwFluid; };

struct InitialState
{
    static constexpr double pressure = BenchmarkCommon::initialPressure;
    static constexpr double temperature = CommonFluid::temperature;
    // Replaces nC10 by the four BSB hydrocarbon lumps while preserving the
    // common benchmark's 20/80 H2O/hydrocarbon overall-composition split.
    inline static constexpr std::array<double, CommonFluid::N> overallComposition{
        0.2000000000000000, 0.0000000100000000,
        0.3340577418658802, 0.2735382352230156,
        0.1333747266780503, 0.0590292862330539};
};

struct Dissolution
{
    static constexpr int component = CommonFluid::co2Component;
    static constexpr double waterMolarMass = 0.018015;
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
    static constexpr bool enableSnesStagnationGuard = true;
    static constexpr int snesStagnationMinimumIterations = 8;
    static constexpr int snesStagnationWindow = 5;
    static constexpr double snesStagnationRelativeImprovement = 1.0e-4;
};

struct Time
{
    static constexpr int numberOfSteps = 50;
    static constexpr double dtDays = 36.525;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = 1.0e-8;
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
    inline static constexpr const char *name = LabScw3D::name;
    using Model = LabScw3D::Model;
    using Grid = LabScw3D::Grid;
    using Fluid = LabScw3D::PrFluid;
    using InitialState = LabScw3D::InitialState;
    using Dissolution = LabScw3D::Dissolution;
    using Land = LabScw3D::Land;
    using Adsorption = LabScw3D::Adsorption;
    using Numerics = LabScw3D::Numerics;
    using Time = LabScw3D::Time;
    using Output = LabScw3D::Output;
    using Rock = LabScw3D::Rock;
};

} // namespace LabScw3D
