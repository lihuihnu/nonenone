/**
 * @file ma2021_case_test.cpp
 * @brief 单元测试：验证 `ma2021_case` 的核心语义、边界条件和回归行为。
 */
#include "../../../case/ma2021_5c_three_phase_reservoir/case_config.hpp"
#include "../../../case/ma2021_5c_three_phase_reservoir/well_config.hpp"

#include <case/well_factory.hpp>
#include <case/case_validation.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/physics/accumulation.hpp>
#include <natural/physics/face_flux.hpp>
#include <natural/physics/well_source.hpp>
#include <natural/state/cell_property_evaluator.hpp>
#include <natural/state/state_codec.hpp>
#include <natural/state/three_phase_equilibrium.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>

namespace
{

using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells,
    CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption,
    CaseConfig::Model::enableLandTrapping,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;
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
        throw std::runtime_error(
            message + ": actual=" + std::to_string(actual) +
            ", expected=" + std::to_string(expected));
}

void relativeNear(double actual, double expected, double tolerance, const std::string &message)
{
    const double scale = std::max(1.0e-12, std::abs(expected));
    if (std::abs(actual - expected) > tolerance * scale)
        throw std::runtime_error(
            message + ": actual=" + std::to_string(actual) +
            ", expected=" + std::to_string(expected));
}

void checkExactLiteratureInputs()
{
    static_assert(Indices::fullyCompositionalThreePhase,
                  "Ma case must use the full O/G/W formulation");
    static_assert(Indices::numComponents == 5,
                  "Ma Table-7 benchmark is a five-component mixture");

    MPMC::cases::validateCaseConfig<Indices, CaseConfig::Config>();

    near(CaseConfig::InitialState::pressure / CaseConfig::bar, 13.79, 1.0e-13,
         "Ma Table-9 pressure");
    near(CaseConfig::InitialState::temperature, 366.5, 1.0e-13,
         "Ma Table-9 temperature");

    const Composition expectedZ{0.10, 0.10, 0.20, 0.40, 0.20};
    for (std::size_t i = 0; i < expectedZ.size(); ++i)
        near(CaseConfig::InitialState::overallComposition[i], expectedZ[i], 1.0e-13,
             "Ma Table-7 overall composition");

    const Composition expectedTc{647.3, 190.6, 507.5, 622.1, 718.6};
    const Composition expectedPcBar{220.47, 46.00, 32.89, 25.34, 18.49};
    const Composition expectedOmega{0.344, 0.008, 0.275, 0.444, 0.651};
    const Composition expectedMw{0.018, 0.016, 0.086, 0.134, 0.206};
    for (std::size_t i = 0; i < expectedZ.size(); ++i)
    {
        near(CaseConfig::Fluid::criticalTemperature[i], expectedTc[i], 1.0e-13,
             "Ma Table-7 Tc");
        near(CaseConfig::Fluid::criticalPressure[i] / CaseConfig::bar,
             expectedPcBar[i], 1.0e-13, "Ma Table-7 Pc");
        near(CaseConfig::Fluid::acentricFactor[i], expectedOmega[i], 1.0e-13,
             "Ma Table-7 acentric factor");
        near(CaseConfig::Fluid::molarMass[i], expectedMw[i], 1.0e-13,
             "Ma Table-7 molecular weight");
    }

    // Selected Table-8 entries protect the published non-aqueous BIP matrix.
    near(CaseConfig::Fluid::binaryInteraction[0][1], 0.4850, 1.0e-13,
         "Ma Table-8 H2O-C1 BIP");
    near(CaseConfig::Fluid::binaryInteraction[0][4], 0.4800, 1.0e-13,
         "Ma Table-8 H2O-C15 BIP");
    near(CaseConfig::Fluid::binaryInteraction[2][3], 0.002866, 1.0e-13,
         "Ma Table-8 C6-C10 BIP");
    near(CaseConfig::Fluid::binaryInteraction[2][4], 0.010970, 1.0e-13,
         "Ma Table-8 C6-C15 BIP");
    near(CaseConfig::Fluid::binaryInteraction[3][4], 0.002657, 1.0e-13,
         "Ma Table-8 C10-C15 BIP");

    near(CaseConfig::LiteratureReference::composition[0][0], 0.00622, 1.0e-13,
         "Ma Table-9 H2O in oil-rich phase");
    near(CaseConfig::LiteratureReference::composition[1][0], 0.059951, 1.0e-13,
         "Ma Table-9 H2O in gas-rich phase");
    near(CaseConfig::LiteratureReference::composition[2][1], 0.000004, 1.0e-13,
         "Ma Table-9 C1 in water-rich phase");
    near(CaseConfig::LiteratureReference::phaseMoleFraction[0], 0.73102, 1.0e-13,
         "Ma Table-9 oil phase fraction");
    near(CaseConfig::LiteratureReference::phaseMoleFraction[1], 0.07813, 1.0e-13,
         "Ma Table-9 gas phase fraction");
    near(CaseConfig::LiteratureReference::phaseMoleFraction[2], 0.19085, 1.0e-13,
         "Ma Table-9 water phase fraction");
}

void checkOrdinaryPrThreePhaseComparison()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::Config>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    const auto flash = equilibrium.flashPTZ(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);

    require(flash.converged, "Ma P-T-z ordinary-PR comparison must converge");
    require(flash.presence.bits() == MPMC::PhasePresence::allBits,
            "Ma P-T-z comparison must contain oil/gas/water-rich phases");

    near(std::accumulate(flash.phaseMoleFraction.begin(), flash.phaseMoleFraction.end(), 0.0),
         1.0, 1.0e-10, "phase mole fractions must close");
    near(std::accumulate(flash.saturation.begin(), flash.saturation.end(), 0.0),
         1.0, 1.0e-10, "phase saturations must close");

    Composition reconstructed{};
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        require(flash.phaseMoleFraction[p] > 0.0 && flash.saturation[p] > 0.0,
                "all three Ma comparison phases must have positive amounts");
        for (int component = 0; component < Indices::numComponents; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            require(std::isfinite(flash.composition[p][c]) && flash.composition[p][c] >= 0.0,
                    "phase composition must remain finite and non-negative");
            reconstructed[c] += flash.phaseMoleFraction[p] * flash.composition[p][c];
        }
    }
    for (std::size_t c = 0; c < reconstructed.size(); ++c)
        near(reconstructed[c], CaseConfig::InitialState::overallComposition[c], 2.0e-8,
             "Ma comparison P-T-z material balance");

    // True mutual solubility: water enters both hydrocarbon phases, and methane
    // enters the aqueous phase.  This is the defining requirement of the new
    // literature case and explicitly excludes the removed Yang formulation.
    require(flash.composition[0][CaseConfig::Fluid::waterComponent] > 1.0e-3,
            "oil-rich phase must contain finite H2O");
    require(flash.composition[1][CaseConfig::Fluid::waterComponent] > 1.0e-2,
            "gas-rich phase must contain finite H2O");
    require(flash.composition[2][CaseConfig::Fluid::methaneComponent] > 1.0e-7,
            "water-rich phase must contain finite methane");
    require(flash.composition[2][CaseConfig::Fluid::waterComponent] > 0.99,
            "third liquid must be water-rich");

    // The ordinary-PR comparison should retain the correct magnitude of the
    // two most visible mutual-solubility markers from Table 9.  We deliberately
    // do NOT assert the phase fractions because Ma et al. use aqueous-specific
    // BIPs and a modified H2O alpha function that v23 does not yet implement.
    relativeNear(flash.composition[0][CaseConfig::Fluid::waterComponent],
                 CaseConfig::LiteratureReference::composition[0][0], 0.10,
                 "ordinary PR H2O-in-oil should remain close to Ma Table 9");
    relativeNear(flash.composition[1][CaseConfig::Fluid::waterComponent],
                 CaseConfig::LiteratureReference::composition[1][0], 0.10,
                 "ordinary PR H2O-in-gas should remain close to Ma Table 9");
}

void checkReservoirKernelPath()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::Config>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    const auto flash = equilibrium.flashPTZ(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);

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
                "each Ma conserved component needs finite positive accumulation");
        componentMass += value;
    }
    double phaseMass = 0.0;
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        require(properties.density[p] > 0.0 && std::isfinite(properties.density[p]),
                "Ma phase density must be finite and positive");
        require(properties.viscosity[p] > 0.0 && std::isfinite(properties.viscosity[p]),
                "Ma phase viscosity must be finite and positive");
        require(properties.mobility[p] > 0.0 && std::isfinite(properties.mobility[p]),
                "Ma phase mobility must be finite and positive");
        phaseMass += properties.porosity * properties.density[p] * properties.saturation[p];
    }
    near(componentMass, phaseMass, 2.0e-10,
         "Ma component accumulation must equal total three-phase fluid mass");
    near(accumulation.waterMass, 0.0, 1.0e-14,
         "full Ma model must not create duplicate independent-water accumulation");

    auto exteriorPrimary = primary;
    exteriorPrimary[Indices::Primary::pressure] =
        CaseConfig::InitialState::pressure - 0.05 * CaseConfig::bar;
    const auto exteriorState =
        MPMC::CellStateCodec<Indices>::decode(exteriorPrimary, phaseState);
    const auto exteriorProperties =
        evaluator.evaluate(exteriorState, CaseConfig::Rock::porosity);
    const auto flux = MPMC::computeFaceMassFlux<Indices>(
        state, properties, exteriorState, exteriorProperties,
        1.0e-12, 0.0, 1.0, 1.0, 0, 1);
    for (double value : flux.component)
        require(std::isfinite(value), "Ma five-component face flux must remain finite");
    require(std::accumulate(flux.component.begin(), flux.component.end(), 0.0) > 0.0,
            "pressure gradient must create a nonzero Ma component flux");
    near(flux.water, 0.0, 1.0e-14,
         "full Ma model must not create a separate water face flux");

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
    std::array<double, Indices::numPhases> injectionPhaseFraction{};
    injectionPhaseFraction[static_cast<std::size_t>(Indices::Phase::vapor)] = 1.0;
    std::array<double, Indices::numComponents> injectionMassFraction{};
    injectionMassFraction[static_cast<std::size_t>(CaseConfig::Fluid::methaneComponent)] = 1.0;

    const auto injector = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Injector, 1.0e-12, 14.0 * CaseConfig::bar,
        phasePressure, density, mobility, fluid.surfaceDensity,
        injectionPhaseFraction, injectionMassFraction, properties.massFraction);
    require(injector.componentMassSource[CaseConfig::Fluid::methaneComponent] > 0.0,
            "Ma wrapper C1 injector must add methane mass");
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        if (component == CaseConfig::Fluid::methaneComponent)
            continue;
        near(injector.componentMassSource[static_cast<std::size_t>(component)],
             0.0, 1.0e-14, "pure-C1 injector must not directly add another component");
    }

    const auto producer = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Producer, 1.0e-12, 13.50 * CaseConfig::bar,
        phasePressure, density, mobility, fluid.surfaceDensity,
        {}, {}, properties.massFraction);
    for (double value : producer.componentMassSource)
        require(value < 0.0 && std::isfinite(value),
                "Ma BHP producer must remove every initially present component");

    near(WellConfig::wells[0].target * CaseConfig::secondsPerDay, 1.0, 1.0e-13,
         "Ma wrapper injector target must equal 1 std m3/day");
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

        checkExactLiteratureInputs();
        checkOrdinaryPrThreePhaseComparison();
        checkReservoirKernelPath();
        std::cout << "Ma 2021 five-component true O/G/W mutual-solubility preflight: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Ma 2021 five-component true O/G/W mutual-solubility preflight: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
