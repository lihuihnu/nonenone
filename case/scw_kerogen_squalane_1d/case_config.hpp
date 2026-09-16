#pragma once

#include "../scw_kerogen_common/benchmark_common.hpp"
#include "../scw_kerogen_common/calibrated_binary_properties.hpp"

#include <indices/model_config.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

/**
 * @file case_config.hpp
 * @brief 第一阶段全组分等温超临界水驱单重质拟组分配置。
 */
namespace CaseConfig
{

inline constexpr char name[] = "scw_kerogen_squalane_1d";

struct Model
{
    static constexpr int numberOfComponents = 2;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase;
};

using Grid = ScwKerogen1D::Grid;
using Rock = ScwKerogen1D::Rock;
using Numerics = ScwKerogen1D::Numerics;
using Time = ScwKerogen1D::Time;
using Output = ScwKerogen1D::Output;

struct Fluid
{
    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int waterComponent = 0;
    static constexpr int heavyComponent = 1;
    inline static constexpr std::array<const char *, N> componentNames{
        "H2O", "Heavy_squalane"};

    // Squalane is branched C30H62.  Tc/Pc/omega follow the existing
    // Teratani-based strict-SCW input; kij is the separate target-window
    // coexistence-composition fit recorded in calibrated_binary_properties.hpp.
    inline static constexpr std::array<double, N> criticalTemperature{
        647.096, 855.9};
    inline static constexpr std::array<double, N> criticalPressure{
        22.064e6, 0.7164e6};
    // PR-derived placeholders used by the existing flash tool.  Vc is unused
    // by unshifted PR fugacity but is used by LBC viscosity, so absolute heavy-
    // phase viscosity remains a screening prediction until Vc/viscosity is fit.
    inline static constexpr std::array<double, N> criticalVolume{
        7.49587808839913e-5, 3.05355324646748e-3};
    inline static constexpr std::array<double, N> acentricFactor{
        0.3443, 1.255};
    inline static constexpr std::array<double, N> molarMass{
        0.01801528, 0.4228};
    inline static constexpr std::array<std::array<double, N>, N>
        binaryInteraction{{
            {{0.0, 0.053233659543451335}},
            {{0.053233659543451335, 0.0}}
        }};

    static double cubicBinaryInteractionCoefficient(
        int i, int j, double temperatureK)
    {
        return i == j ? 0.0
                      : ScwKerogenCalibration::prSqualaneKij(temperatureK);
    }

    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;
    // PR76 is retained to match the existing H2O-squalane benchmark.
    static constexpr int eosModelFlag = 1;
    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr double eosU = 2.4142135623730951;
    static constexpr double eosW = -0.4142135623730951;
    static constexpr double temperature = ScwKerogen1D::temperature;

    inline static constexpr std::array<double, 3> surfaceDensity{
        800.0, 120.0, 360.0};
    inline static constexpr std::array<double, 3> viscosity{
        1.0e-3, 5.0e-5, 7.0e-5};
    static constexpr double waterViscosity = 7.0e-5;
    static constexpr double waterFormationVolumeFactor = 1.0;
    static constexpr bool useIapws2008AqueousViscosity = true;
    static constexpr double aqueousViscosityMaximumSoluteMoleFraction = 0.02;
};

struct InitialState
{
    static constexpr double pressure = ScwKerogen1D::initialPressure;
    static constexpr double temperature = ScwKerogen1D::temperature;
    // The selected z_H2O=0.60 lies on the one-phase side of the measured
    // hydrocarbon-rich endpoint x_H2O=0.911 at 653.2 K and 27.74 MPa.  The
    // experiment-supported baseline therefore starts as one hydrocarbon-rich
    // phase and allows the injected H2O front to create a water-rich phase.
    inline static constexpr std::array<double, Fluid::N> overallComposition{
        ScwKerogenCalibration::SqualaneTargetState::initialOverallH2OMoleFraction,
        1.0 - ScwKerogenCalibration::SqualaneTargetState::initialOverallH2OMoleFraction};
};

struct Dissolution
{
    static constexpr int component = Fluid::heavyComponent;
    static constexpr double waterMolarMass = 0.01801528;
    static constexpr double salinityMolality = 0.0;
    static constexpr double initialWaterCO2MoleFraction = 0.0;
};
struct Land { static constexpr double constant = 0.0; };
struct Adsorption
{
    static constexpr double rockDensity = 2650.0;
    inline static constexpr std::array<double, Fluid::N> thetaMax{};
    inline static constexpr std::array<double, Fluid::N> coefficient{};
    static constexpr double standardPressure = 101325.0;
    static constexpr double standardTemperature = 288.15;
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
