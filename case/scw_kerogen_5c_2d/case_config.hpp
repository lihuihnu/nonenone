#pragma once
#include "../scw_kerogen_common/benchmark_common.hpp"
#include <indices/model_config.hpp>
#include <natural/thermo/thermodynamic_model.hpp>
#include <array>

namespace ScwKerogen5C2D {

struct Model {
    static constexpr int numberOfComponents = 5;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase;
};

struct Grid {
    inline static constexpr const char *meshDirectory = "";
    static constexpr int nx = 60, ny = 20, nz = 1;
    static constexpr double lx = 0.30, ly = 0.10, lz = 0.010;
};

struct Rock {
    static constexpr double kx = 1e-12, ky = 1e-12, kz = 1e-12;
    static constexpr double porosity = 0.25;
};

struct Fluid {
    static constexpr std::size_t N = 5;
    static constexpr int waterComponent = 0;
    static constexpr int gasolineComponent = 1;
    static constexpr int dieselComponent = 2;
    static constexpr int middleComponent = 3;
    static constexpr int heavyComponent = 4;

    inline static constexpr std::array<int, 2>
        producerLighteningLightComponents{1, 2};
    inline static constexpr std::array<int, 1>
        producerLighteningHeavyComponents{4};

    inline static constexpr std::array<const char *, N> componentNames{
        "H2O", "OIL_GASOLINE", "OIL_DIESEL", "OIL_MIDDLE", "OIL_HEAVY"};

    inline static constexpr std::array<double, N> criticalTemperature{
        647.096, 515.231, 736.423, 870.437, 982.864};
    inline static constexpr std::array<double, N> criticalPressure{
        22.064e6, 3.192179e6, 1.751155e6, 1.138075e6, 0.808714e6};
    inline static constexpr std::array<double, N> criticalVolume{
        55.948074534e-6, 358.490e-6, 824.073e-6,
        1318.853e-6, 1793.890e-6};
    inline static constexpr std::array<double, N> acentricFactor{
        0.3443, 0.269328, 0.583473, 0.896554, 1.215676};
    inline static constexpr std::array<double, N> molarMass{
        0.01801528, 0.082392, 0.215583, 0.387349, 0.660132};

    // 653.15 K PR screening priors from
    // scw_kerogen_lumped_flow_study/pr_parameters/pr_binary_matrix_screening.csv.
    // They are intentionally retained as conditional screening inputs and are
    // NOT promoted to experimentally validated production BIPs.
    inline static constexpr std::array<std::array<double, N>, N>
        binaryInteraction{{
            {{0.0,          0.5000000000, 0.6662345314, 0.2398345519, 0.2398345519}},
            {{0.5000000000, 0.0,          0.0,          0.0,          0.0}},
            {{0.6662345314, 0.0,          0.0,          0.0,          0.0}},
            {{0.2398345519, 0.0,          0.0,          0.0,          0.0}},
            {{0.2398345519, 0.0,          0.0,          0.0,          0.0}}
        }};

    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;
    static constexpr int eosModelFlag = 1;
    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr double eosU = 2.4142135623730951;
    static constexpr double eosW = -0.4142135623730951;
    static constexpr double temperature = 653.15;

    // Legacy phase-property defaults are only fallback/reference values.
    // Production compositional properties use the common LBC/IAPWS paths.
    inline static constexpr std::array<double, 3> surfaceDensity{
        800.0, 120.0, 360.0};
    inline static constexpr std::array<double, 3> viscosity{
        8.0e-4, 4.0e-5, 5.0e-5};
    static constexpr double waterViscosity = 5.0e-5;
    static constexpr double waterFormationVolumeFactor = 1.0;
    static constexpr bool useIapws2008AqueousViscosity = true;
    static constexpr double aqueousViscosityMaximumSoluteMoleFraction = 0.02;
};

struct InitialState {
    static constexpr double pressure = 28e6;
    static constexpr double temperature = 653.15;
    // Formal five-component BASE anchor: z_H2O=0.20, with the characterized
    // recovered-oil mole ratio preserved in the remaining 0.80.
    inline static constexpr std::array<double, Fluid::N> overallComposition{
        0.20,
        0.0290640,
        0.3254160,
        0.2603368,
        0.1851832
    };
};

struct Dissolution {
    static constexpr int component = Fluid::gasolineComponent;
    static constexpr double waterMolarMass = 0.01801528;
    static constexpr double salinityMolality = 0.0;
    static constexpr double initialWaterCO2MoleFraction = 0.0;
};

struct Land { static constexpr double constant = 0.0; };

struct Adsorption {
    static constexpr double rockDensity = 2650.0;
    static constexpr double standardPressure = 101325.0;
    static constexpr double standardTemperature = 288.15;
    inline static constexpr std::array<double, Fluid::N> thetaMax{};
    inline static constexpr std::array<double, Fluid::N> coefficient{};
};

struct Numerics : ScwKerogen1D::Numerics {
    static constexpr double massResidualScale = 2e-5;
    static constexpr double rateWellResidualFloor = 2.5e-8;
};

struct Time : ScwKerogen1D::Time {
    static constexpr int numberOfSteps = 75; // 4500 s safety horizon
    static constexpr double dtDays = 60.0 / 86400.0;
    static constexpr double maximumDtDays = 2.0 / 86400.0;
    static constexpr double minimumDtDays = 1e-5 / 86400.0;
    static constexpr int maximumRetries = 16;
    static constexpr double targetPVI = 1.0;
};

struct Output : ScwKerogen1D::Output {
    inline static constexpr const char *directory = "./results";
    static constexpr std::size_t every = 25; // field snapshots at 0.5 PVI spacing
    // Numerical screening keeps the H02 effective PV so actual-PVI accounting
    // and flow scale remain directly comparable to the accepted binary study.
    static constexpr double producerEffectivePoreVolumeM3 = 7.5e-5;
};

struct Config {
    inline static constexpr const char *name = "scw_kerogen_5c_2d";
    using Model = ScwKerogen5C2D::Model;
    using Grid = ScwKerogen5C2D::Grid;
    using Rock = ScwKerogen5C2D::Rock;
    using Fluid = ScwKerogen5C2D::Fluid;
    using InitialState = ScwKerogen5C2D::InitialState;
    using Dissolution = ScwKerogen5C2D::Dissolution;
    using Land = ScwKerogen5C2D::Land;
    using Adsorption = ScwKerogen5C2D::Adsorption;
    using Numerics = ScwKerogen5C2D::Numerics;
    using Time = ScwKerogen5C2D::Time;
    using Output = ScwKerogen5C2D::Output;
};

} // namespace ScwKerogen5C2D
