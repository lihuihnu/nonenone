/**
 * @file three_phase_reservoir_case_test.cpp
 * @brief 单元测试：验证 `three_phase_reservoir_case` 的核心语义、边界条件和回归行为。
 */
#include "../../../case/3p4c_pr_reservoir/case_config.hpp"
#include "../../../case/3p4c_pr_reservoir/well_config.hpp"

#include <case/well_factory.hpp>
#include <case/case_validation.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/fluid_system.hpp>
#include <natural/assembly/cell_residual.hpp>
#include <natural/petsc/property_conversion.hpp>
#include <natural/physics/accumulation.hpp>
#include <natural/physics/face_flux.hpp>
#include <natural/physics/well_source.hpp>
#include <natural/state/cell_property_evaluator.hpp>
#include <natural/state/state_codec.hpp>
#include <natural/state/three_phase_equilibrium.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <well/peaceman.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells,
    false,
    false,
    false,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;
using AdIndices = MPMC::ADIndices<ModelConfig>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Composition = std::array<double, Indices::numComponents>;

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string &message)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (std::abs(actual - expected) > tolerance * scale)
    {
        throw std::runtime_error(
            message + ": actual=" + std::to_string(actual) +
            ", expected=" + std::to_string(expected));
    }
}

Eos makeEos()
{
    MPMC::CompositionalMixture<Indices> mixture(
        CaseConfig::Fluid::criticalTemperature,
        CaseConfig::Fluid::criticalPressure,
        CaseConfig::Fluid::criticalVolume,
        CaseConfig::Fluid::acentricFactor,
        CaseConfig::Fluid::molarMass,
        CaseConfig::Fluid::binaryInteraction);

    return Eos(
        CaseConfig::Fluid::eosOmegaA,
        CaseConfig::Fluid::eosOmegaB,
        std::move(mixture),
        CaseConfig::Fluid::eosModelFlag,
        CaseConfig::Fluid::eosU,
        CaseConfig::Fluid::eosW,
        1.0e-30);
}

MPMC::FluidSystem<Indices> makeFluid()
{
    return MPMC::cases::makeFluidSystem<Indices, CaseConfig::Config>();
}

void checkCaseFlash()
{
    MPMC::cases::validateCaseConfig<Indices, CaseConfig::Config>();
    const auto eos = makeEos();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = CaseConfig::Fluid::waterComponent;
    MPMC::CubicThreePhaseFlash<Indices> flash(eos, options);

    const auto result = flash.flash(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);

    require(result.converged, "reservoir-case P-T-z flash must converge");
    require(result.presence.bits() == MPMC::PhasePresence::allBits,
            "reservoir-case initial state must contain O+G+W");

    const std::array<double, 3> betaReference{
        0.165091541711314,
        0.187058670847636,
        0.647849787441050};
    const std::array<double, 3> saturationReference{
        0.377870658521436,
        0.317428566143412,
        0.304700775335152};

    for (std::size_t phase = 0; phase < 3; ++phase)
    {
        near(result.phaseMoleFraction[phase], betaReference[phase], 2.0e-8,
             "case phase mole fraction regression");
        near(result.saturation[phase], saturationReference[phase], 2.0e-8,
             "case saturation regression");
        require(result.saturation[phase] > 0.25,
                "initial O/G/W saturations should all be substantial for the reservoir smoke case");
    }

    // Physical phase identity checks for this selected P-T-z state.
    const auto &oil = result.composition[0];
    const auto &gas = result.composition[1];
    const auto &water = result.composition[2];
    require(oil[CaseConfig::Fluid::heavyComponent] > 0.10,
            "oil-rich phase must contain substantial nC16");
    require(gas[CaseConfig::Fluid::co2Component] > 0.85,
            "gas phase must be CO2-rich");
    require(water[CaseConfig::Fluid::waterComponent] > 0.99,
            "water-rich phase must be H2O-dominated");
    require(oil[CaseConfig::Fluid::waterComponent] > 1.0e-3,
            "oil-rich phase must contain finite H2O");
    require(gas[CaseConfig::Fluid::waterComponent] > 1.0e-3,
            "gas phase must contain finite H2O");
    require(water[CaseConfig::Fluid::co2Component] > 1.0e-4,
            "water-rich phase must contain finite CO2");

    // Reconstruct z from beta*x; this is the defining P-T-z material balance.
    Composition reconstructed{};
    for (int phase = 0; phase < 3; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        for (int component = 0; component < Indices::numComponents; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            reconstructed[c] += result.phaseMoleFraction[p] * result.composition[p][c];
        }
    }
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        near(reconstructed[c], CaseConfig::InitialState::overallComposition[c], 2.0e-9,
             "reservoir-case flash material balance");
    }

    // Fugacity equality for components not pinned at the numerical composition
    // floor in the water-rich phase.
    const auto oilResult = eos.phaseResult(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        oil,
        true);
    const auto gasResult = eos.phaseResult(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        gas,
        false);
    const auto waterResult = eos.phaseResult(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        water,
        true);

    for (int component = 0; component < Indices::numComponents; ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        const double logOg = std::abs(std::log(
            oilResult.fugacity[c] / gasResult.fugacity[c]));
        require(logOg < 2.0e-8, "oil-gas fugacity equality failed");

        if (water[c] > 1.0e-20)
        {
            const double logOw = std::abs(std::log(
                oilResult.fugacity[c] / waterResult.fugacity[c]));
            require(logOw < 2.0e-8, "oil-water fugacity equality failed");
        }
    }
}

void checkCaseFeedsReservoirPhysics()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    const auto flash = equilibrium.flashPTZ(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);
    require(flash.converged && flash.presence.count() == 3,
            "reservoir physics preflight requires the three-phase initial state");

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = CaseConfig::InitialState::pressure;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, flash);

    const auto state = MPMC::CellStateCodec<Indices>::decode(primary, phaseState);
    MPMC::CellPropertyEvaluator<Indices> evaluator(fluid);
    const auto properties = evaluator.evaluate(state, CaseConfig::Rock::porosity);
    const auto accumulation = MPMC::computeFluidAccumulation<Indices>(properties);

    double componentMass = 0.0;
    for (double value : accumulation.componentMass)
    {
        require(value > 0.0 && std::isfinite(value),
                "all four components must have finite positive initial accumulation");
        componentMass += value;
    }

    double phaseMass = 0.0;
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        require(properties.density[p] > 0.0 && std::isfinite(properties.density[p]),
                "all initial phase densities must be finite and positive");
        require(properties.viscosity[p] > 0.0 && std::isfinite(properties.viscosity[p]),
                "all initial phase viscosities must be finite and positive");
        require(properties.mobility[p] > 0.0 && std::isfinite(properties.mobility[p]),
                "all initial phase mobilities must be finite and positive");
        phaseMass += properties.porosity * properties.density[p] * properties.saturation[p];
    }
    near(componentMass, phaseMass, 2.0e-10,
         "component accumulation must equal total three-phase fluid mass");

    // A one-bar pressure drop must produce a finite three-phase component flux.
    // This does not replace PETSc integration, but it verifies that the actual
    // case fluid can enter the same Darcy face-flux kernel used by the reservoir.
    auto exteriorPrimary = primary;
    exteriorPrimary[Indices::Primary::pressure] =
        CaseConfig::InitialState::pressure - CaseConfig::bar;
    const auto exteriorState =
        MPMC::CellStateCodec<Indices>::decode(exteriorPrimary, phaseState);
    const auto exteriorProperties =
        evaluator.evaluate(exteriorState, CaseConfig::Rock::porosity);
    const auto flux = MPMC::computeFaceMassFlux<Indices>(
        state,
        properties,
        exteriorState,
        exteriorProperties,
        1.0e-12,
        0.0,
        1.0,
        1.0,
        0,
        1);
    double totalAbsoluteFlux = 0.0;
    for (double value : flux.component)
    {
        require(std::isfinite(value), "case component face flux must be finite");
        totalAbsoluteFlux += std::abs(value);
    }
    require(totalAbsoluteFlux > 0.0,
            "case pressure gradient must create a nonzero component face flux");
    near(flux.water, 0.0, 1.0e-14,
         "full model must not create a duplicate independent-water face flux");

    // Exercise the actual case injection/production semantics in the common
    // compositional well-source kernel.
    std::array<double, Indices::numPhases> phasePressure{};
    std::array<double, Indices::numPhases> phaseDensity{};
    std::array<double, Indices::numPhases> phaseMobility{};
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        phasePressure[p] = CaseConfig::InitialState::pressure;
        phaseDensity[p] = properties.density[p];
        phaseMobility[p] = properties.mobility[p];
    }

    std::array<double, Indices::numPhases> injectionPhaseFraction{};
    injectionPhaseFraction[static_cast<std::size_t>(Indices::Phase::vapor)] = 1.0;
    std::array<double, Indices::numComponents> injectionMassFraction{};
    injectionMassFraction[static_cast<std::size_t>(CaseConfig::Fluid::co2Component)] = 1.0;

    const auto injector = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Injector,
        1.0e-12,
        205.0 * CaseConfig::bar,
        phasePressure,
        phaseDensity,
        phaseMobility,
        fluid.surfaceDensity,
        injectionPhaseFraction,
        injectionMassFraction,
        properties.massFraction);
    require(injector.componentMassSource[
                static_cast<std::size_t>(CaseConfig::Fluid::co2Component)] > 0.0,
            "pure-CO2 injector must add CO2 mass");
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        if (component == CaseConfig::Fluid::co2Component)
            continue;
        near(injector.componentMassSource[static_cast<std::size_t>(component)],
             0.0, 1.0e-14,
             "pure-CO2 injector must not directly add another component");
    }

    const auto producer = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Producer,
        1.0e-12,
        195.0 * CaseConfig::bar,
        phasePressure,
        phaseDensity,
        phaseMobility,
        fluid.surfaceDensity,
        {},
        {},
        properties.massFraction);
    for (double value : producer.componentMassSource)
        require(value < 0.0 && std::isfinite(value),
                "BHP producer must remove each initially present conserved component");
    near(producer.waterMassSource, 0.0, 1.0e-14,
         "full model must not create a duplicate independent-water well source");

    // Sanity-check that the configured rate is deliberately mild relative to
    // reservoir pore volume; this is a smoke case, not a shock test.
    const double bulkVolume =
        CaseConfig::Grid::lx * CaseConfig::Grid::ly * CaseConfig::Grid::lz;
    const double poreVolume = bulkVolume * CaseConfig::Rock::porosity;
    const double injectedVolumePerDay =
        WellConfig::wells[0].target * CaseConfig::secondsPerDay;
    near(injectedVolumePerDay, 100.0, 1.0e-12,
         "well-config rate conversion must equal 100 m3/day");
    require(injectedVolumePerDay / poreVolume < 1.0e-3,
            "default CO2 rate should be below 0.1% PV/day for the initial smoke test");
}

void checkInitialReservoirResidualScale()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    const auto flash = equilibrium.flashPTZ(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = CaseConfig::InitialState::pressure;
    primary[Indices::Primary::wellPressure] = CaseConfig::InitialState::pressure;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, flash);

    const auto state = MPMC::CellStateCodec<Indices>::decode(primary, phaseState);
    MPMC::CellPropertyEvaluator<Indices> evaluator(fluid);
    const auto properties = evaluator.evaluate(state, CaseConfig::Rock::porosity);

    std::array<double, Indices::numPhases> phasePressure{};
    std::array<double, Indices::numPhases> density{};
    std::array<double, Indices::numPhases> mobility{};
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        phasePressure[p] = CaseConfig::InitialState::pressure;
        density[p] = properties.density[p];
        mobility[p] = properties.mobility[p];
    }

    const double dx = CaseConfig::Grid::lx / CaseConfig::Grid::nx;
    const double dy = CaseConfig::Grid::ly / CaseConfig::Grid::ny;
    const double dz = CaseConfig::Grid::lz / CaseConfig::Grid::nz;
    const double wellIndex = MPMC::verticalPeacemanWellIndex(
        {dx, dy, dz, CaseConfig::Rock::kx, CaseConfig::Rock::ky, 0.10, 0.0});

    std::array<double, Indices::numPhases> injectionPhaseFraction{};
    injectionPhaseFraction[static_cast<std::size_t>(Indices::Phase::vapor)] = 1.0;
    std::array<double, Indices::numComponents> injectionMassFraction{};
    injectionMassFraction[static_cast<std::size_t>(CaseConfig::Fluid::co2Component)] = 1.0;

    const auto injector = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Injector, wellIndex,
        WellConfig::wells[0].initialBhp, phasePressure, density, mobility,
        fluid.surfaceDensity, injectionPhaseFraction, injectionMassFraction,
        properties.massFraction);
    const auto producer = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Producer, wellIndex,
        WellConfig::wells[1].initialBhp, phasePressure, density, mobility,
        fluid.surfaceDensity, {}, {}, properties.massFraction);

    MPMC::FaceMassFlux<Indices, double> zeroFlux{};
    const double cellVolume = dx * dy * dz;
    const double dt = CaseConfig::Time::dtDays * CaseConfig::secondsPerDay;
    const auto injectorResidual = MPMC::assembleCellResidual<Indices>(
        state, properties, properties, zeroFlux,
        injector.componentMassSource, injector.waterMassSource,
        cellVolume, dt, 1.0);
    const auto producerResidual = MPMC::assembleCellResidual<Indices>(
        state, properties, properties, zeroFlux,
        producer.componentMassSource, producer.waterMassSource,
        cellVolume, dt, 1.0);

    const double injectorRateResidual =
        MPMC::selectControlledRate<Indices>(
            MPMC::WellControl::TotalRate, injector.surfacePhaseRate) -
        WellConfig::wells[0].target;

    double normSquared = injectorRateResidual * injectorRateResidual;
    for (int equation = 0; equation < Indices::numEquations; ++equation)
    {
        if (equation == Indices::Equation::wellControl)
            continue;
        const double injectorValue =
            injectorResidual.value[static_cast<std::size_t>(equation)];
        const double producerValue =
            producerResidual.value[static_cast<std::size_t>(equation)];
        require(std::isfinite(injectorValue) && std::isfinite(producerValue),
                "initial production-case well residuals must remain finite");
        normSquared += injectorValue * injectorValue + producerValue * producerValue;
    }

    const double twoWellResidualNorm = std::sqrt(normSquared);
    require(twoWellResidualNorm > 0.1 && twoWellResidualNorm < 10.0,
            "initial 3p4c residual should be controlled by physical well sources, not a thermodynamic spike");
}

void checkInitialTraceComponentResidualRegression()
{
    using AdValue = typename AdIndices::ValueType;

    auto fluid = MPMC::cases::makeFluidSystem<AdIndices, CaseConfig::Config>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<AdIndices> equilibrium(fluid);
    const auto flash = equilibrium.flashPTZ(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);
    require(flash.converged && flash.presence.count() == 3,
            "trace-component regression requires the production O+G+W initial flash");

    std::array<double, AdIndices::numPrimaryVariables> primary{};
    primary[AdIndices::Primary::pressure] = CaseConfig::InitialState::pressure;
    primary[AdIndices::Primary::wellPressure] = CaseConfig::InitialState::pressure;
    MPMC::PhaseStateData<AdIndices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, flash);

    const auto state = MPMC::CellStateCodec<AdIndices>::decode(primary, phaseState);
    const std::size_t heavy = static_cast<std::size_t>(CaseConfig::Fluid::heavyComponent);
    require(MPMC::scalarValue(state.aqueousMoleFraction[heavy]) >= 0.0,
            "dependent trace component must not become negative after primary decoding");
    require(MPMC::scalarValue(state.aqueousMoleFraction[heavy]) <=
                MPMC::NaturalNumerics::phaseEquilibriumTraceComposition,
            "water-rich nC16 must remain on the numerical trace-composition boundary");

    MPMC::CellPropertyEvaluator<AdIndices> evaluator(fluid);
    const auto properties = evaluator.evaluate(state, AdValue(CaseConfig::Rock::porosity));
    const auto previous = MPMC::scalarizeCellProperties<AdIndices>(properties);
    MPMC::FaceMassFlux<AdIndices, AdValue> zeroFlux{};
    std::array<AdValue, AdIndices::numComponents> zeroWell{};

    const double cellVolume =
        (CaseConfig::Grid::lx / CaseConfig::Grid::nx) *
        (CaseConfig::Grid::ly / CaseConfig::Grid::ny) *
        (CaseConfig::Grid::lz / CaseConfig::Grid::nz);
    const auto residual = MPMC::assembleCellResidual<AdIndices>(
        state, properties, previous, zeroFlux, zeroWell, AdValue(0.0),
        cellVolume, CaseConfig::Time::dtDays * CaseConfig::secondsPerDay, 1.0);

    double maximumResidual = 0.0;
    double maximumJacobianEntry = 0.0;
    for (const auto &equation : residual.value)
    {
        const double value = MPMC::scalarValue(equation);
        require(std::isfinite(value),
                "production initial local residual must remain finite");
        maximumResidual = std::max(maximumResidual, std::abs(value));
        for (int derivative = 0; derivative < AdIndices::numPrimaryVariables; ++derivative)
        {
            const double jacobian = equation.derivative(derivative);
            require(std::isfinite(jacobian),
                    "production initial local Jacobian must remain finite");
            maximumJacobianEntry = std::max(maximumJacobianEntry, std::abs(jacobian));
        }
    }

    require(maximumResidual < 1.0e-6,
            "equilibrated production initial state must not contain a large local residual");
    require(maximumJacobianEntry < 1.0e12,
            "trace-component handling must not create an extreme local Jacobian entry");

    const int waterHeavyEquation =
        AdIndices::Equation::waterFugacity[heavy];
    near(MPMC::scalarValue(residual.value[static_cast<std::size_t>(waterHeavyEquation)]),
         0.0, 1.0e-14,
         "trace water-rich nC16 row must be a zero-composition boundary equation");
}

} // namespace

int main()
{
    try
    {
        {
            std::vector<MPMC::WellPerforation<int>> perforations{{0, 1.0e-12}};
            const auto well = MPMC::cases::makeWell<Indices, int>(
                WellConfig::wells.front(), 0, std::move(perforations));
            require(well.id == WellConfig::wells.front().id,
                    "case well factory must preserve the configured well id");
            require(well.control == MPMC::cases::detail::toNaturalControl(
                        WellConfig::wells.front().control),
                    "case well factory must preserve the configured well control");
        }

        checkCaseFlash();
        checkCaseFeedsReservoirPhysics();
        checkInitialReservoirResidualScale();
        checkInitialTraceComponentResidualRegression();
        std::cout << "H2O-CO2-CH4-nC16 reservoir case preflight: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "H2O-CO2-CH4-nC16 reservoir case preflight: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
