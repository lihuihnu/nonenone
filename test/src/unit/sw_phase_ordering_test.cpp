/**
 * @file sw_phase_ordering_test.cpp
 * @brief 回归二维算例中 SW 油气相槽位反转问题。
 */
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/fluid_system.hpp>
#include <natural/state/state_codec.hpp>
#include <natural/state/three_phase_equilibrium.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace {
using Config = MPMC::CompositionalModelConfig<
    3, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Equilibrium = MPMC::FullyCompositionalThreePhaseEquilibrium<Indices>;
using Composition = std::array<double, 3>;

constexpr std::array<double, 3> criticalTemperature{647.096, 304.1282, 617.70};
constexpr std::array<double, 3> criticalPressure{22.064e6, 7.3773e6, 2.103e6};
constexpr std::array<double, 3> criticalVolume{5.5948e-5, 9.4118e-5, 6.10e-4};
constexpr std::array<double, 3> acentricFactor{0.3443, 0.22394, 0.4920};
constexpr std::array<double, 3> molarMass{0.018015268, 0.0440098, 0.14228168};
constexpr std::array<std::array<double, 3>, 3> binaryInteraction{{
    {{0.0, 0.1896, 0.5000}},
    {{0.1896, 0.0, 0.1141}},
    {{0.5000, 0.1141, 0.0}}
}};

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

Eos makeSw()
{
    MPMC::CompositionalMixture<Indices> mixture(
        criticalTemperature, criticalPressure, criticalVolume,
        acentricFactor, molarMass, binaryInteraction);
    Eos eos(
        0.45724, 0.07780, std::move(mixture), 1,
        2.4142135623730951, -0.4142135623730951, 1.0e-30);

    Eos::SoreideWhitsonOptions sw;
    sw.waterComponent = 0;
    sw.salinityMolality = 0.0;
    sw.aqueousWaterBip[0] = [](double, double) { return 0.0; };
    sw.aqueousWaterBip[1] = [](double, double) { return -0.06609; };
    sw.aqueousWaterBip[2] = [](double, double) { return -0.14229; };
    eos.configureSoreideWhitson(std::move(sw));
    return eos;
}

MPMC::FluidSystem<Indices> makeFluid()
{
    MPMC::FluidSystem<Indices> fluid(
        {730.0, 1.8, 985.4040020947351},
        {2.4e-4, 1.86e-5, 4.69091e-4},
        makeSw(),
        {"H2O", "CO2", "nC10"},
        {"Oil", "Gas", "Water"},
        333.15);
    fluid.configureFullyCompositionalThreePhase(0);
    return fluid;
}

void writeComposition(
    std::array<double, Indices::numPrimaryVariables>& primary,
    const std::array<int, Indices::numIndependentCompositionsPerPhase>& indices,
    const Composition& composition)
{
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        primary[static_cast<std::size_t>(
            indices[static_cast<std::size_t>(component)])] =
            composition[static_cast<std::size_t>(component)];
    }
}

void checkBenchmarkFrontState()
{
    // Convex combination of the three coexistence compositions recorded at
    // cell 490 of the original 0.1-PVI SW run.  Before the fix, the returned
    // public oil slot is CO2-rich and the gas slot is nC10-rich.
    constexpr Composition co2Rich{
        0.006082530464619983, 0.9923058453430816, 0.0016116241922984};
    constexpr Composition nc10Rich{
        0.002326074097273210, 0.3900437088477006, 0.6076302170550262};
    constexpr Composition waterRich{
        0.9900629231333399, 0.009937076866523781, 1.36318761812e-13};
    constexpr std::array<double, 3> beta{0.20, 0.50, 0.30};
    Composition overall{};
    for (std::size_t component = 0; component < overall.size(); ++component)
        overall[component] =
            beta[0] * co2Rich[component] +
            beta[1] * nc10Rich[component] +
            beta[2] * waterRich[component];

    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    const auto eos = makeSw();
    const Flash flash(eos, options);
    const auto result = flash.flash(4.813700167512818e6, 333.15, overall);

    require(result.converged, "SW benchmark-front flash must converge");
    require(result.presence.bits() == MPMC::PhasePresence::allBits,
            "SW benchmark-front flash must contain O+G+W");
    require(result.composition[2][0] > 0.99,
            "SW benchmark-front water slot must remain water-rich");
    require(result.composition[1][1] > result.composition[0][1],
            "SW benchmark-front gas must be richer in CO2 than oil");
    require(result.composition[0][2] > result.composition[1][2],
            "SW benchmark-front oil must be richer in nC10 than gas");

    // Removing the physical gas produces a metastable O+W restricted
    // equilibrium.  The stability test must still recover the true missing
    // gas; the reference-equilibrium certification must not suppress a real
    // phase transition.
    MPMC::PhasePresence oilWater(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit);
    const auto reduced = flash.flashRestricted(
        4.813700167512818e6, 333.15, overall, oilWater, result.composition);
    require(reduced.converged && reduced.presence.bits() == oilWater.bits(),
            "SW benchmark-front restricted O+W state must converge");
    const auto stability = flash.stabilityTest(
        4.813700167512818e6, 333.15, overall,
        reduced.presence, reduced.composition);
    require(stability.valid && stability.missingPhaseUnstable[1],
            "SW benchmark-front TPD must recover the true missing gas");
    require(stability.trialSum[1] > 1.1,
            "SW benchmark-front missing gas must have a strong TPD signal");
    require(stability.incipientComposition[1][1] > 0.9,
            "SW benchmark-front incipient gas must be CO2-rich, not aqueous");
}

void checkNewtonStateCanonicalization()
{
    const auto fluid = makeFluid();
    const Equilibrium equilibrium(fluid);
    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 4.813700167512818e6;
    primary[Indices::Primary::liquidSaturation] = 0.07266240554205435;
    primary[Indices::Primary::vaporSaturation] = 0.7226364550775557;
    primary[Indices::Primary::waterSaturation] = 0.20470113938038997;

    const Composition co2Rich{
        0.006082530464619983, 0.9923058453430816, 0.0016116241922984};
    const Composition nc10Rich{
        0.002326074097273210, 0.3900437088477006, 0.6076302170550262};
    const Composition waterRich{
        0.9900629231333399, 0.009937076866523781, 1.36318761812e-13};
    writeComposition(primary, Indices::Primary::liquidComposition, co2Rich);
    writeComposition(primary, Indices::Primary::vaporComposition, nc10Rich);
    writeComposition(primary, Indices::Primary::waterComposition, waterRich);

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence::all();
    equilibrium.updateSecondary(primary, phaseState);
    const Composition overallBefore = phaseState.overallComposition;
    equilibrium.updatePhaseState(primary, phaseState);

    require(
        primary[Indices::Primary::vaporComposition[1]] >
            primary[Indices::Primary::liquidComposition[1]],
        "SW Newton state must map the CO2-rich phase to gas");
    require(
        primary[Indices::Primary::liquidComposition[0]] <
            primary[Indices::Primary::vaporComposition[0]],
        "SW Newton state must map the lower-water nonaqueous phase to oil");
    for (std::size_t component = 0; component < overallBefore.size(); ++component)
    {
        require(
            std::abs(phaseState.overallComposition[component] -
                     overallBefore[component]) < 2.0e-13,
            "SW Newton phase relabeling must preserve overall composition");
    }
}

void checkTransientNewtonStabilityUsesRestrictedEquilibriumReference()
{
    const auto fluid = makeFluid();
    const Equilibrium equilibrium(fluid);
    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 5.2e6;
    primary[Indices::Primary::liquidSaturation] = 0.8;
    primary[Indices::Primary::vaporSaturation] = 0.0;
    primary[Indices::Primary::waterSaturation] = 0.2;

    // A reproducible OW Newton intermediate state.  Its phase records are
    // deliberately not yet at interphase fugacity equilibrium, but their
    // beta-weighted overall composition has a stable physical O+W flash.
    // The old SW TPD path used the raw oil record as the reference chemical
    // potential and reported a false missing Gas with trialSum=1.0012749;
    // that "incipient gas" was 99.9962 mol% H2O (the aqueous basin).
    constexpr Composition oilIterate{
        0.001780516781684561,
        0.034561874244895971,
        0.96365760897341946};
    constexpr Composition waterIterate{
        0.99904170781500368,
        0.00095829218399018942,
        1.0061779184449603e-12};
    constexpr double betaOil = 0.26023139184001487;
    constexpr double betaWater = 0.73976860815998513;
    Composition overall{};
    for (std::size_t component = 0; component < overall.size(); ++component)
        overall[component] =
            betaOil * oilIterate[component] +
            betaWater * waterIterate[component];

    writeComposition(primary, Indices::Primary::liquidComposition, oilIterate);
    writeComposition(primary, Indices::Primary::vaporComposition, overall);
    writeComposition(primary, Indices::Primary::waterComposition, waterIterate);

    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    const Flash flash(fluid.eos, options);
    const std::array<Composition, 3> transientCompositions{
        oilIterate, overall, waterIterate};
    const MPMC::PhasePresence oilWater(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit);
    const auto stability = flash.stabilityTest(
        5.2e6, 333.15, overall, oilWater, transientCompositions);
    require(stability.valid,
            "SW transient OW stability certification must be valid");
    require(!stability.missingPhaseUnstable[1],
            "SW transient OW TPD must reject the aqueous false-gas basin");

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = oilWater;
    phaseState.phaseMoleFraction = {betaOil, 0.0, betaWater};
    phaseState.overallComposition = overall;

    const auto result = equilibrium.updatePhaseState(primary, phaseState);
    require(!result.missingPhaseUnstable,
            "SW transient OW iterate must not create a false missing gas");
    require(phaseState.phasePresence.bits() ==
                (MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit),
            "SW transient OW iterate must remain O+W");

    // No phase transition was certified, so the Newton unknowns for the active
    // phases must remain untouched; only secondary state may be refreshed.
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        require(std::abs(
                    primary[static_cast<std::size_t>(
                        Indices::Primary::liquidComposition[c])] -
                    oilIterate[c]) < 1.0e-14,
                "SW stability check must not overwrite oil Newton iterate");
        require(std::abs(
                    primary[static_cast<std::size_t>(
                        Indices::Primary::waterComposition[c])] -
                    waterIterate[c]) < 1.0e-14,
                "SW stability check must not overwrite water Newton iterate");
    }
}

void checkDependentCompositionCancellationBoundary()
{
    using AdIndices = MPMC::ADIndices<Config>;
    using Codec = MPMC::CellStateCodec<AdIndices>;

    typename Codec::PrimaryArray primary{};
    primary[AdIndices::Primary::pressure] = 5.16e6;
    primary[AdIndices::Primary::liquidSaturation] = 0.8;
    primary[AdIndices::Primary::vaporSaturation] = 0.0;
    primary[AdIndices::Primary::waterSaturation] = 0.2;

    constexpr double stalledWaterNc10 = 9.7810648469476291e-14;
    primary[static_cast<std::size_t>(AdIndices::Primary::waterComposition[0])] =
        0.999;
    primary[static_cast<std::size_t>(AdIndices::Primary::waterComposition[1])] =
        0.001 - stalledWaterNc10;

    MPMC::PhaseStateData<AdIndices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit);

    const auto boundaryState = Codec::decode(primary, phaseState);
    const auto &dependent = boundaryState.aqueousMoleFraction[2];
    require(std::abs(dependent.value()) < 1.0e-30,
            "dependent SW trace composition must snap to zero boundary");
    require(dependent.derivative(AdIndices::Primary::waterComposition[0]) == -1.0,
            "dependent boundary must retain first composition derivative");
    require(dependent.derivative(AdIndices::Primary::waterComposition[1]) == -1.0,
            "dependent boundary must retain second composition derivative");

    // The cancellation treatment is not a global 1e-13 trace threshold: a
    // clearly resolvable dependent value above the measured boundary remains a
    // normal positive composition and therefore continues to use thermodynamic
    // equilibrium equations.
    constexpr double resolvableWaterNc10 = 2.0e-13;
    primary[static_cast<std::size_t>(AdIndices::Primary::waterComposition[1])] =
        0.001 - resolvableWaterNc10;
    const auto resolvedState = Codec::decode(primary, phaseState);
    require(resolvedState.aqueousMoleFraction[2].value() >
                MPMC::NaturalNumerics::dependentCompositionCancellationBoundary,
            "resolvable dependent composition must not be snapped to zero");
}

} // namespace

int main()
{
    try {
        checkBenchmarkFrontState();
        checkNewtonStateCanonicalization();
        checkTransientNewtonStabilityUsesRestrictedEquilibriumReference();
        checkDependentCompositionCancellationBoundary();
        std::cout << "SW benchmark-front phase ordering: ALL PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SW benchmark-front phase ordering: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
