/**
 * @file panfili2025_case3_test.cpp
 * @brief 单元测试：验证 `panfili2025_case3` 的核心语义、边界条件和回归行为。
 */
#include "../../../case/panfili2025_case3_fullphysics/case_config.hpp"
#include "../../../case/panfili2025_case3_fullphysics/case_fluid.hpp"
#include "../../../case/panfili2025_case3_fullphysics/well_config.hpp"

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

void checkPaperInputsAndDeclaredAdaptations()
{
    static_assert(Indices::fullyCompositionalThreePhase,
                  "Panfili 2025 case #3 requires the live-water O/G/W formulation");
    static_assert(Indices::numComponents == 10,
                  "Panfili case #3 is nine non-water components plus H2O");
    static_assert(CaseConfig::Fluid::thermodynamicModel ==
                      MPMC::CubicThermodynamicModel::SoreideWhitson,
                  "Panfili case #3 must use Soreide-Whitson thermodynamics");

    MPMC::cases::validateCaseConfig<Indices, CaseConfig::Config>();

    require(CaseConfig::Grid::nx == 47 && CaseConfig::Grid::ny == 99 && CaseConfig::Grid::nz == 15,
            "paper grid dimensions must be retained");
    near(CaseConfig::Grid::dx / CaseConfig::foot, 328.0, 1.0e-13,
         "paper average dx");
    near(CaseConfig::Grid::dy / CaseConfig::foot, 328.0, 1.0e-13,
         "paper average dy");
    near(CaseConfig::Grid::dz / CaseConfig::foot, 8.0, 1.0e-13,
         "paper average dz");
    near(CaseConfig::Rock::kx / CaseConfig::mD, 815.0, 1.0e-13,
         "paper average horizontal permeability");
    near(CaseConfig::Rock::porosity, 0.15, 1.0e-13,
         "paper average porosity");
    near(CaseConfig::LiteratureReference::initialPressure / CaseConfig::psi,
         3000.0, 1.0e-13, "paper initial pressure");
    near(CaseConfig::LiteratureReference::referenceDepth / CaseConfig::foot,
         2910.0, 1.0e-13, "paper gas-water-contact reference depth");
    near((CaseConfig::Fluid::temperature - 273.15) * 9.0 / 5.0 + 32.0,
         200.0, 1.0e-12, "paper reservoir temperature");
    near(CaseConfig::LiteratureReference::connateWaterSaturation,
         0.30, 1.0e-13, "Section 7.3 connate-water saturation");

    const std::array<double, 9> expectedZ{
        0.01210, 0.01940, 0.65990, 0.08690, 0.05910,
        0.12970, 0.02745, 0.00515, 0.00030};
    const std::array<double, 9> expectedMw{
        44.01, 28.013, 16.043, 30.07, 44.097,
        66.869, 107.779, 198.562, 335.198};
    const std::array<double, 9> expectedTcR{
        548.46, 227.16, 343.08, 549.77, 665.64,
        806.54, 838.11, 1058.04, 1291.89};
    const std::array<double, 9> expectedPcPsi{
        1071.33, 492.31, 667.78, 708.34, 618.70,
        514.93, 410.75, 247.56, 160.42};
    const std::array<double, 9> expectedOmega{
        0.22500, 0.04000, 0.01300, 0.09860, 0.15240,
        0.21575, 0.31230, 0.55670, 0.91692};

    near(std::accumulate(expectedZ.begin(), expectedZ.end(), 0.0), 1.0, 1.0e-13,
         "Panfili Table-3 non-water feed must close");
    for (std::size_t j = 0; j < expectedZ.size(); ++j)
    {
        near(CaseConfig::LiteratureReference::hydrocarbonComposition[j], expectedZ[j],
             1.0e-13, "Panfili Table-3 mole fraction");
        const std::size_t c = j + 1;
        near(CaseConfig::Fluid::molarMass[c] * 1000.0, expectedMw[j], 2.0e-13,
             "Panfili Table-3 molecular weight");
        near(CaseConfig::Fluid::criticalTemperature[c] * 9.0 / 5.0, expectedTcR[j],
             2.0e-13, "Panfili Table-3 critical temperature");
        near(CaseConfig::Fluid::criticalPressure[c] / CaseConfig::psi, expectedPcPsi[j],
             2.0e-13, "Panfili Table-3 critical pressure");
        near(CaseConfig::Fluid::acentricFactor[c], expectedOmega[j], 1.0e-13,
             "Panfili Table-3 acentric factor");
    }

    // Protect the explicitly documented SPE3 closure used only because the
    // Panfili paper does not print its hydrocarbon-hydrocarbon BIP matrix.
    near(CaseConfig::Fluid::binaryInteraction[1][2], -0.0200, 1.0e-13,
         "SPE3 closure CO2-N2 BIP");
    near(CaseConfig::Fluid::binaryInteraction[1][3], 0.1000, 1.0e-13,
         "SPE3 closure CO2-C1 BIP");
    near(CaseConfig::Fluid::binaryInteraction[3][6], 0.09281, 1.0e-13,
         "SPE3 closure C1-C4-6 BIP");
    near(CaseConfig::Fluid::binaryInteraction[3][9], 0.13920, 1.0e-13,
         "SPE3 closure C1-C7+3 BIP");
}

void checkFullPhysicsInitializationAndPhaseAppearance()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::Config>();
    require(fluid.eos.usesSoreideWhitson(), "production fluid factory must select SW");

    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    const auto initial = equilibrium.flashPTZ(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);

    require(initial.converged, "Panfili initial SW flash must converge");
    require(initial.presence.bits() ==
                (MPMC::PhasePresence::gasBit | MPMC::PhasePresence::waterBit),
            "3000-psi initial state must be undersaturated G+W");
    near(initial.saturation[Indices::Phase::liquid], 0.0, 1.0e-12,
         "initial oil saturation");
    near(initial.saturation[Indices::Phase::vapor], 0.70, 2.0e-9,
         "derived initial gas saturation");
    near(initial.saturation[Indices::Phase::water], 0.30, 2.0e-9,
         "paper-target initial connate water saturation");

    Composition reconstructed{};
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        for (int component = 0; component < Indices::numComponents; ++component)
            reconstructed[static_cast<std::size_t>(component)] +=
                initial.phaseMoleFraction[p] *
                initial.composition[p][static_cast<std::size_t>(component)];
    }
    for (std::size_t c = 0; c < reconstructed.size(); ++c)
        near(reconstructed[c], CaseConfig::InitialState::overallComposition[c], 3.0e-8,
             "initial SW P-T-z material balance");

    const double nonWaterTotal = 1.0 - CaseConfig::InitialState::overallComposition[0];
    for (std::size_t j = 0; j < CaseConfig::LiteratureReference::hydrocarbonComposition.size(); ++j)
        near(CaseConfig::InitialState::overallComposition[j + 1] / nonWaterTotal,
             CaseConfig::LiteratureReference::hydrocarbonComposition[j], 3.0e-12,
             "derived 30%-Sw initialization must preserve Panfili Table-3 non-water ratios");

    // Full live-water physics: water vaporization and non-water dissolution
    // are both present, rather than using the restricted Section-7.1 switches.
    require(initial.composition[Indices::Phase::vapor][CaseConfig::Fluid::waterComponent] > 1.0e-3,
            "full physics must vaporize water into the gas phase");
    require(initial.composition[Indices::Phase::water][CaseConfig::Fluid::co2Component] > 1.0e-5,
            "full physics must dissolve CO2 into the aqueous phase");
    require(initial.composition[Indices::Phase::water][CaseConfig::Fluid::methaneComponent] > 1.0e-4,
            "full physics must permit hydrocarbon dissolution into water");

    const auto at2000 = equilibrium.flashPTZ(
        2000.0 * CaseConfig::psi,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);
    require(at2000.converged &&
                at2000.presence.bits() ==
                    (MPMC::PhasePresence::gasBit | MPMC::PhasePresence::waterBit),
            "fluid should remain above the condensate boundary near 2000 psi");

    const auto at1800 = equilibrium.flashPTZ(
        1800.0 * CaseConfig::psi,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);
    require(at1800.converged, "1800-psi SW flash must converge");
    require(at1800.presence.bits() == MPMC::PhasePresence::allBits,
            "depletion must be capable of producing the intended G+W -> O+G+W transition");
    require(at1800.saturation[Indices::Phase::liquid] > 1.0e-4,
            "new condensate phase must have positive oil saturation");
}

void checkScheduleAndPureCo2Injection()
{
    require(WellConfig::wells.size() == 2,
            "paper schedule must keep exactly the same two physical crest wells");
    require(WellConfig::wells[0].completion.i != WellConfig::wells[1].completion.i ||
                WellConfig::wells[0].completion.j != WellConfig::wells[1].completion.j,
            "the two physical wells require distinct BHP representative cells");

    near(WellConfig::producerTarget(),
         WellConfig::mscfPerDayToM3PerSecond(50000.0), 1.0e-13,
         "paper depletion rate per well");
    near(CaseConfig::LiteratureReference::producerMinimumBhp / CaseConfig::psi,
         435.0, 1.0e-13,
         "paper minimum producer BHP");
    near(WellConfig::injectorTarget(),
         WellConfig::mscfPerDayToM3PerSecond(90000.0), 1.0e-13,
         "paper injection rate per well");
    near(CaseConfig::LiteratureReference::injectorMaximumBhp / CaseConfig::psi,
         2850.0, 1.0e-13,
         "paper maximum injector BHP");

    const double y31 = WellConfig::years(31.0);
    const double y32 = WellConfig::years(32.0);
    const double y33 = WellConfig::years(33.0);
    const double y34 = WellConfig::years(34.0);
    const double y35 = WellConfig::years(35.0);
    const double y49 = WellConfig::years(49.0);
    const double y50 = WellConfig::years(50.0);
    const double y51 = WellConfig::years(51.0);

    require(WellConfig::stageAt(y31) == WellConfig::OperatingStage::Depletion,
            "producer must be active before year 32");
    require(WellConfig::stageAt(y32) == WellConfig::OperatingStage::Depletion,
            "31->32 y backward-Euler step must remain a production step");
    require(WellConfig::stageAt(y33) == WellConfig::OperatingStage::Idle,
            "year 33 must be inside the two-year idle period");
    require(WellConfig::stageAt(y34) == WellConfig::OperatingStage::Idle,
            "year 34 is the schedule boundary; injection begins on the next step");
    require(WellConfig::stageAt(y35) == WellConfig::OperatingStage::Injection,
            "34->35 y step must use pure-CO2 injection");
    require(WellConfig::stageAt(y49) == WellConfig::OperatingStage::Injection,
            "injector must be active in the Fig.25 injection window");
    require(WellConfig::stageAt(y50) == WellConfig::OperatingStage::Injection,
            "49->50 y step must remain an injection step");
    require(WellConfig::stageAt(y51) == WellConfig::OperatingStage::Closed,
            "comparison wrapper must close the wells after the ~year-50 injection period");
}

void checkReservoirKernelPath()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::Config>();
    Panfili2025CaseFluid::applyRockFluidApproximation(fluid);
    near(fluid.waterRelativePermeability(0.30), 0.0, 1.0e-14,
         "paper-target connate water must be immobile in the wrapper closure");
    require(fluid.waterRelativePermeability(0.60) > 0.0,
            "water relative permeability must activate above connate water");
    require(fluid.gasRelativePermeability(0.35) > 0.0,
            "gas relative permeability must remain active");
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
        require(std::isfinite(value) && value >= 0.0,
                "Panfili component accumulation must remain finite/non-negative");
        componentMass += value;
    }
    require(componentMass > 0.0, "Panfili total fluid accumulation must be positive");
    near(accumulation.waterMass, 0.0, 1.0e-14,
         "live-water formulation must not duplicate an independent water equation");

    double phaseMass = 0.0;
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        if (properties.saturation[p] <= 0.0)
            continue;
        require(properties.density[p] > 0.0 && std::isfinite(properties.density[p]),
                "active Panfili phase density must be finite/positive");
        require(properties.viscosity[p] > 0.0 && std::isfinite(properties.viscosity[p]),
                "active Panfili phase viscosity must be finite/positive");
        phaseMass += properties.porosity * properties.density[p] * properties.saturation[p];
    }
    near(componentMass, phaseMass, 5.0e-10,
         "Panfili component accumulation must equal total active-phase mass");

    auto exteriorPrimary = primary;
    exteriorPrimary[Indices::Primary::pressure] -= 5.0 * CaseConfig::psi;
    const auto exteriorState = MPMC::CellStateCodec<Indices>::decode(exteriorPrimary, phaseState);
    const auto exteriorProperties = evaluator.evaluate(exteriorState, CaseConfig::Rock::porosity);
    const auto flux = MPMC::computeFaceMassFlux<Indices>(
        state, properties, exteriorState, exteriorProperties,
        1.0e-12, 0.0, 1.0, 1.0, 0, 1);
    for (double value : flux.component)
        require(std::isfinite(value), "Panfili 10-component face flux must remain finite");
    near(flux.water, 0.0, 1.0e-14,
         "live-water model must not create a separate water face flux");

    std::array<double, Indices::numPhases> phasePressure{};
    std::array<double, Indices::numPhases> density{};
    std::array<double, Indices::numPhases> mobility{};
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        phasePressure[p] = CaseConfig::InitialState::pressure;
        density[p] = std::max(properties.density[p], 1.0);
        mobility[p] = std::max(properties.mobility[p], 0.0);
    }
    std::array<double, Indices::numPhases> injectionPhaseFraction{};
    injectionPhaseFraction[static_cast<std::size_t>(Indices::Phase::vapor)] = 1.0;
    std::array<double, Indices::numComponents> injectionMassFraction{};
    injectionMassFraction[CaseConfig::Fluid::co2Component] = 1.0;

    const auto injector = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Injector, 1.0e-12, 3100.0 * CaseConfig::psi,
        phasePressure, density, mobility, fluid.surfaceDensity,
        injectionPhaseFraction, injectionMassFraction, properties.massFraction);
    require(injector.componentMassSource[CaseConfig::Fluid::co2Component] > 0.0,
            "pure CO2 injector must add CO2 mass");
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        if (component == CaseConfig::Fluid::co2Component)
            continue;
        near(injector.componentMassSource[static_cast<std::size_t>(component)],
             0.0, 1.0e-14, "pure CO2 injector must not directly add another component");
    }
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

        checkPaperInputsAndDeclaredAdaptations();
        checkFullPhysicsInitializationAndPhaseAppearance();
        checkScheduleAndPureCo2Injection();
        checkReservoirKernelPath();
        std::cout << "Panfili 2025 case #3 full-physics live-water preflight: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Panfili 2025 case #3 full-physics live-water preflight: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
