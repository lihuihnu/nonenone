/**
 * @file three_phase_flash_test.cpp
 * @brief 单元测试：验证 `three_phase_flash` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <natural/fluid_system.hpp>
#include <natural/state/three_phase_equilibrium.hpp>
#include <natural/state/state_codec.hpp>
#include <natural/state/cell_property_evaluator.hpp>
#include <natural/physics/accumulation.hpp>
#include <natural/physics/face_flux.hpp>
#include <natural/physics/well_source.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace
{

using Config = MPMC::CompositionalModelConfig<
    4,
    true,
    false,
    false,
    false,
    false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Composition = std::array<double, Indices::numComponents>;
using Eos = MPMC::CubicEquationOfState<Indices>;

static_assert(std::is_same_v<
    MPMC::CubicThreePhaseFlash<Indices>,
    MPMC::PengRobinsonThreePhaseFlash<Indices>>,
    "Legacy three-phase flash alias must remain source compatible.");

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
    // Synthetic CO2/CH4/n-C16/H2O system.  The deliberately large water-
    // hydrocarbon BICs create a water-rich liquid while still allowing finite
    // mutual solubility in all three PR phases.
    MPMC::CompositionalMixture<Indices> mixture(
        {304.2, 190.6, 717.0, 647.3},
        {73.8e5, 46.0e5, 14.2e5, 220.5e5},
        {9.4e-5, 9.9e-5, 9.0e-4, 5.6e-5},
        {0.225, 0.008, 0.742, 0.344},
        {0.0440098, 0.016043, 0.22644, 0.01801528},
        {{0.0, 0.1000, 0.1250, 0.1896},
         {0.1000, 0.0, 0.0780, 0.4850},
         {0.1250, 0.0780, 0.0, 0.5000},
         {0.1896, 0.4850, 0.5000, 0.0}});

    // Ordinary Peng-Robinson constants requested for the first implementation.
    return Eos(
        0.4572355,
        0.0779691,
        std::move(mixture),
        1,
        2.414213562373095,
        -0.414213562373095,
        1.0e-30);
}

MPMC::FluidSystem<Indices> makeFluid()
{
    auto fluid = MPMC::FluidSystem<Indices>(
        {800.0, 20.0, 1000.0},
        {1.0e-3, 1.0e-5, 1.0e-4},
        makeEos(),
        {"CO2", "CH4", "nC16", "H2O"},
        {"Oil", "Gas", "Water"},
        350.0);
    fluid.configureFullyCompositionalThreePhase(3);
    fluid.gasRelativePermeability = [](double s) { return s * s; };
    fluid.waterRelativePermeability = [](double s) { return s * s; };
    fluid.threePhaseOilRelativePermeability = [](double, double so, double) { return so * so; };
    return fluid;
}

void checkComposition(const Composition &x, const std::string &name)
{
    double sum = 0.0;
    for (double value : x)
    {
        require(std::isfinite(value) && value >= 0.0,
                name + " contains invalid mole fraction");
        sum += value;
    }
    near(sum, 1.0, 1.0e-10, name + " mole fractions must sum to one");
}

void testThreePhasePTzFlash()
{
    auto eos = makeEos();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 3;
    MPMC::CubicThreePhaseFlash<Indices> flash(eos, options);

    constexpr double pressure = 50.0e5;
    constexpr double temperature = 350.0;
    const Composition z{0.75, 0.025, 0.025, 0.20};

    const auto result = flash.flash(pressure, temperature, z);
    require(result.converged, "PR three-phase PTz flash must converge");
    require(result.presence.bits() == MPMC::PhasePresence::allBits,
            "known PTz state must contain oil/gas/water-rich phases");

    near(std::accumulate(result.phaseMoleFraction.begin(), result.phaseMoleFraction.end(), 0.0),
         1.0, 1.0e-10, "phase mole fractions must close");
    near(std::accumulate(result.saturation.begin(), result.saturation.end(), 0.0),
         1.0, 1.0e-10, "phase saturations must close");

    for (int p = 0; p < 3; ++p)
    {
        require(result.phaseMoleFraction[static_cast<std::size_t>(p)] > 0.0,
                "three-phase mole fraction must be positive");
        require(result.saturation[static_cast<std::size_t>(p)] > 0.0,
                "three-phase saturation must be positive");
        require(result.compressibility[static_cast<std::size_t>(p)] > 0.0,
                "three-phase compressibility must be positive");
        checkComposition(result.composition[static_cast<std::size_t>(p)],
                         "phase composition");
    }

    // Reconstruct the user-supplied overall composition from beta*x.
    Composition reconstructed{};
    for (int p = 0; p < 3; ++p)
        for (int c = 0; c < Indices::numComponents; ++c)
            reconstructed[static_cast<std::size_t>(c)] +=
                result.phaseMoleFraction[static_cast<std::size_t>(p)] *
                result.composition[static_cast<std::size_t>(p)][static_cast<std::size_t>(c)];
    for (int c = 0; c < Indices::numComponents; ++c)
        near(reconstructed[static_cast<std::size_t>(c)], z[static_cast<std::size_t>(c)],
             2.0e-8, "PTz material balance");

    // Demonstrate true mutual partitioning: H2O is finite in oil/gas and CO2 is
    // finite in the water-rich phase; no Henry-only immiscible-water shortcut.
    require(result.composition[0][3] > 1.0e-5, "oil-rich phase must contain H2O");
    require(result.composition[1][3] > 1.0e-5, "gas phase must contain H2O");
    require(result.composition[2][0] > 1.0e-6, "water-rich phase must contain CO2");
    require(result.composition[2][3] > 0.9, "third liquid must be water-rich");

    // Fugacity equality is the thermodynamic closure for every component.
    const auto oil = eos.phaseResult(pressure, temperature, result.composition[0], true);
    const auto gas = eos.phaseResult(pressure, temperature, result.composition[1], false);
    const auto water = eos.phaseResult(pressure, temperature, result.composition[2], true);
    near(eos.compressibility(pressure, temperature, result.composition[0], MPMC::CompositionalPhase::Oil),
         oil.compressibility, 1.0e-14, "PR compressibility-only oil root");
    near(eos.compressibility(pressure, temperature, result.composition[1], MPMC::CompositionalPhase::Gas),
         gas.compressibility, 1.0e-14, "PR compressibility-only gas root");
    near(eos.compressibility(pressure, temperature, result.composition[2], MPMC::CompositionalPhase::Water),
         water.compressibility, 1.0e-14, "PR compressibility-only water root");
    for (int c = 0; c < Indices::numComponents; ++c)
    {
        const std::size_t i = static_cast<std::size_t>(c);
        const auto relativeError = [](double a, double b) {
            return std::abs(a - b) / std::max({1.0, std::abs(a), std::abs(b)});
        };
        require(relativeError(oil.fugacity[i], gas.fugacity[i]) < 2.0e-6,
                "oil-gas component fugacity mismatch");
        // The n-C16 aqueous mole fraction can sit at the numerical composition
        // floor; skip a meaningless relative check only for floor-bound traces.
        if (result.composition[2][i] > 10.0 * options.compositionFloor)
            require(relativeError(oil.fugacity[i], water.fugacity[i]) < 2.0e-6,
                    "oil-water component fugacity mismatch");
    }
}

void testPhaseRemovalAndStabilityReappearance()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    constexpr double pressure = 50.0e5;
    const Composition z{0.75, 0.025, 0.025, 0.20};
    const auto flash = equilibrium.flashPTZ(pressure, fluid.temperature, z);
    require(flash.converged && flash.presence.count() == 3,
            "phase-switch regression requires a three-phase reference state");

    // Exercise all three disappearance directions.  Each trial starts from the
    // same physical three-phase P-T-z state, drives exactly one active phase
    // saturation negative, and verifies that the reduced phase set is tested
    // for stability against the unchanged overall composition.
    const std::array<int, 3> saturationIndex{
        Indices::Primary::liquidSaturation,
        Indices::Primary::vaporSaturation,
        Indices::Primary::waterSaturation};
    const std::array<const char *, 3> phaseName{
        "oil-rich", "gas", "water-rich"};

    for (std::size_t missing = 0; missing < saturationIndex.size(); ++missing)
    {
        std::array<double, Indices::numPrimaryVariables> primary{};
        primary[Indices::Primary::pressure] = pressure;
        MPMC::PhaseStateData<Indices> phaseState;
        equilibrium.assignFlashResult(primary, phaseState, flash);

        primary[static_cast<std::size_t>(saturationIndex[missing])] = -1.0e-6;
        equilibrium.updatePhaseState(primary, phaseState);

        require(phaseState.phasePresence.count() == 3,
                std::string("stability test must reintroduce unstable missing ") +
                    phaseName[missing] + " phase");
        require(primary[static_cast<std::size_t>(saturationIndex[missing])] > 0.0,
                std::string("reintroduced ") + phaseName[missing] +
                    " phase must have positive saturation");

        for (int c = 0; c < Indices::numComponents; ++c)
            near(phaseState.overallComposition[static_cast<std::size_t>(c)],
                 z[static_cast<std::size_t>(c)], 2.0e-8,
                 "phase switching must preserve overall composition");
    }
}

void testTraceAppearanceHoldPreventsImmediateDeletion()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    const Composition z{0.75, 0.025, 0.025, 0.20};
    const auto flash = equilibrium.flashPTZ(50.0e5, fluid.temperature, z);
    require(flash.converged && flash.presence.count() == 3,
            "trace appearance-hold regression requires a three-phase reference state");

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 50.0e5;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, flash);

    constexpr double traceWaterSaturation = 2.5e-7;
    const double gasSaturation = primary[Indices::Primary::vaporSaturation];
    primary[Indices::Primary::waterSaturation] = traceWaterSaturation;
    primary[Indices::Primary::liquidSaturation] =
        1.0 - gasSaturation - traceWaterSaturation;
    phaseState.phaseSuppression.add(MPMC::CompositionalPhase::Water);

    const auto held =
        equilibrium.updatePhaseState(primary, phaseState);
    require(!held.phaseRemoved,
            "a certified trace reappearance inside the probe band must not be deleted again");
    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Water),
            "trace water phase must remain active while its appearance hold is active");
    require(phaseState.phaseSuppression.contains(MPMC::CompositionalPhase::Water),
            "trace reappearance hold must persist while saturation remains inside the probe band");
    near(primary[Indices::Primary::waterSaturation],
         traceWaterSaturation, 1.0e-12,
         "appearance hold must not inflate the physical trace saturation");

    // Once the phase grows beyond the probe band it becomes an ordinary active
    // phase and the transition-memory bit is cleared.
    constexpr double grownWaterSaturation = 2.0e-4;
    primary[Indices::Primary::waterSaturation] = grownWaterSaturation;
    primary[Indices::Primary::liquidSaturation] =
        1.0 - gasSaturation - grownWaterSaturation;
    const auto grown =
        equilibrium.updatePhaseState(primary, phaseState);
    require(!grown.phaseRemoved,
            "a grown reappeared phase must remain active");
    require(!phaseState.phaseSuppression.contains(MPMC::CompositionalPhase::Water),
            "appearance hold must clear after saturation leaves the probe band");

    // The hold is not a residual-saturation floor: a non-positive Newton
    // saturation can still trigger the normal restricted/stability transition.
    phaseState.phaseSuppression.add(MPMC::CompositionalPhase::Water);
    primary[Indices::Primary::waterSaturation] = -1.0e-6;
    primary[Indices::Primary::liquidSaturation] =
        1.0 - gasSaturation + 1.0e-6;
    const auto negative =
        equilibrium.updatePhaseState(primary, phaseState);
    require(negative.phaseRemoved,
            "appearance hold must never protect a zero/negative phase saturation");
}

void testFlowPhysicsUsesAllThreeCompositionalPhases()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    const Composition z{0.75, 0.025, 0.025, 0.20};
    const auto flash = equilibrium.flashPTZ(50.0e5, fluid.temperature, z);
    require(flash.converged && flash.presence.count() == 3,
            "flow-physics regression requires a three-phase state");

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 50.0e5;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, flash);

    const auto state = MPMC::CellStateCodec<Indices>::decode(primary, phaseState);
    MPMC::CellPropertyEvaluator<Indices> evaluator(fluid);
    const auto properties = evaluator.evaluate(state, 0.25);
    const auto accumulation = MPMC::computeFluidAccumulation<Indices>(properties);

    double totalComponentMass = 0.0;
    for (double value : accumulation.componentMass)
    {
        require(value > 0.0 && std::isfinite(value),
                "every conserved component must have finite positive accumulation");
        totalComponentMass += value;
    }
    double totalPhaseMass = 0.0;
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        totalPhaseMass += properties.porosity * properties.density[p] * properties.saturation[p];
    }
    near(totalComponentMass, totalPhaseMass, 2.0e-10,
         "sum of component accumulation must equal all-phase fluid mass");
    require(accumulation.componentMass[3] > 0.0,
            "H2O must be conserved as a normal component");
    near(accumulation.waterMass, 0.0, 1.0e-14,
         "full model must not create a duplicate independent-water accumulation");

    auto exteriorPrimary = primary;
    exteriorPrimary[Indices::Primary::pressure] = 49.0e5;
    const auto exteriorState = MPMC::CellStateCodec<Indices>::decode(exteriorPrimary, phaseState);
    const auto exteriorProperties = evaluator.evaluate(exteriorState, 0.25);
    const auto flux = MPMC::computeFaceMassFlux<Indices>(
        state, properties, exteriorState, exteriorProperties,
        1.0e-12, 0.0, 1.0, 1.0, 0, 1);
    require(flux.component[3] > 0.0,
            "water-rich phase must transport H2O through component flux");
    near(flux.water, 0.0, 1.0e-14,
         "full model must not create a separate water face flux");

    std::array<double, Indices::numPhases> phasePressure{};
    std::array<double, Indices::numPhases> density{};
    std::array<double, Indices::numPhases> mobility{};
    for (int phase = 0; phase < Indices::numPhases; ++phase)
    {
        const std::size_t p = static_cast<std::size_t>(phase);
        phasePressure[p] = primary[Indices::Primary::pressure];
        density[p] = properties.density[p];
        mobility[p] = properties.mobility[p];
    }
    std::array<double, Indices::numPhases> injectionPhaseFraction{};
    std::array<double, Indices::numComponents> injectionComponentMassFraction{};
    const auto well = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Producer,
        1.0e-12,
        45.0e5,
        phasePressure,
        density,
        mobility,
        fluid.surfaceDensity,
        injectionPhaseFraction,
        injectionComponentMassFraction,
        properties.massFraction);
    require(well.componentMassSource[3] < 0.0,
            "producer must remove H2O through the common component source");
    near(well.waterMassSource, 0.0, 1.0e-14,
         "full model must not create a duplicate independent-water well source");
}

class FailingPhaseUpdateFlash final
{
public:
    using Result = MPMC::ThreePhaseFlashResult<Indices>;
    using StabilityResult = MPMC::ThreePhaseStabilityResult<Indices>;
    using EosType = MPMC::CubicEquationOfState<Indices>;

    FailingPhaseUpdateFlash(
        const EosType &,
        MPMC::ThreePhaseFlashOptions options = {})
        : options_(options)
    {
    }

    [[nodiscard]] const MPMC::ThreePhaseFlashOptions &options() const noexcept
    {
        return options_;
    }

    [[nodiscard]] Result flash(
        double, double, Composition) const
    {
        return {};
    }

    [[nodiscard]] Result flashRestricted(
        double, double, Composition, MPMC::PhasePresence,
        const std::array<Composition, 3> &) const
    {
        return {};
    }

    [[nodiscard]] StabilityResult stabilityTest(
        double, double, const Composition &, MPMC::PhasePresence active,
        const std::array<Composition, 3> &) const
    {
        StabilityResult result;
        result.valid = true;
        result.stable = false;
        result.testedPresence = active;
        return result;
    }

private:
    MPMC::ThreePhaseFlashOptions options_{};
};

class InvalidStabilityPhaseUpdateFlash final
{
public:
    using Result = MPMC::ThreePhaseFlashResult<Indices>;
    using StabilityResult = MPMC::ThreePhaseStabilityResult<Indices>;
    using EosType = MPMC::CubicEquationOfState<Indices>;

    InvalidStabilityPhaseUpdateFlash(
        const EosType &,
        MPMC::ThreePhaseFlashOptions options = {})
        : options_(options)
    {
    }

    [[nodiscard]] const MPMC::ThreePhaseFlashOptions &options() const noexcept
    {
        return options_;
    }

    [[nodiscard]] Result flash(double, double, Composition) const
    {
        return {};
    }

    [[nodiscard]] Result flashRestricted(
        double, double, Composition, MPMC::PhasePresence active,
        const std::array<Composition, 3> &) const
    {
        Result result;
        result.converged = true;
        result.presence = active;
        return result;
    }

    [[nodiscard]] StabilityResult stabilityTest(
        double, double, const Composition &, MPMC::PhasePresence active,
        const std::array<Composition, 3> &) const
    {
        StabilityResult result;
        result.valid = false;
        result.stable = false;
        result.testedPresence = active;
        return result;
    }

private:
    MPMC::ThreePhaseFlashOptions options_{};
};

void testRecoverablePhaseUpdateFailureIsTransactional()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<
        Indices, FailingPhaseUpdateFlash> equilibrium(fluid);

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 50.0e5;
    primary[Indices::Primary::liquidSaturation] = -1.0e-6;
    primary[Indices::Primary::vaporSaturation] = 0.70;
    primary[Indices::Primary::waterSaturation] = 0.30;

    const std::array<Composition, 3> phaseComposition{{
        Composition{0.70, 0.10, 0.05, 0.15},
        Composition{0.20, 0.60, 0.05, 0.15},
        Composition{0.05, 0.05, 0.05, 0.85}}};
    for (MPMC::CompositionalPhase phase : {
             MPMC::CompositionalPhase::Oil,
             MPMC::CompositionalPhase::Gas,
             MPMC::CompositionalPhase::Water})
    {
        MPMC::three_phase_detail::writeCompositionToPrimary<Indices>(
            primary, phaseComposition[static_cast<std::size_t>(MPMC::phaseIndex(phase))], phase);
    }

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence::all();
    phaseState.phaseMoleFraction = {0.25, 0.45, 0.30};
    for (int c = 0; c < Indices::numComponents; ++c)
    {
        const std::size_t i = static_cast<std::size_t>(c);
        phaseState.overallComposition[i] =
            0.25 * phaseComposition[0][i] +
            0.45 * phaseComposition[1][i] +
            0.30 * phaseComposition[2][i];
    }

    const auto originalPrimary = primary;
    const auto originalState = phaseState;
    const auto result = equilibrium.updatePhaseState(primary, phaseState);

    require(result.recoverableFailure(),
            "failed restricted + unrestricted flash must report recoverable failure");
    require(result.phaseRemoved,
            "recoverable failure must record the attempted active-set removal");
    require(result.restrictedFlashFailed && result.unrestrictedFlashFailed,
            "recoverable failure must identify both failed flash paths");
    require(primary == originalPrimary,
            "failed phase update must restore the complete primary state");
    require(phaseState.phasePresence.bits() == originalState.phasePresence.bits(),
            "failed phase update must restore phase presence");
    require(phaseState.phaseSuppression.bits() == originalState.phaseSuppression.bits(),
            "failed phase update must restore phase suppression history");
    require(phaseState.phaseMoleFraction == originalState.phaseMoleFraction,
            "failed phase update must restore phase mole fractions");
    require(phaseState.overallComposition == originalState.overallComposition,
            "failed phase update must restore overall composition");
}


void testInvalidStabilityFailureIsDistinguishedFromInstability()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<
        Indices, InvalidStabilityPhaseUpdateFlash> equilibrium(fluid);

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 50.0e5;
    primary[Indices::Primary::liquidSaturation] = 0.0;
    primary[Indices::Primary::vaporSaturation] = 0.70;
    primary[Indices::Primary::waterSaturation] = 0.30;

    const std::array<Composition, 3> phaseComposition{{
        Composition{0.70, 0.10, 0.05, 0.15},
        Composition{0.20, 0.60, 0.05, 0.15},
        Composition{0.05, 0.05, 0.05, 0.85}}};
    for (MPMC::CompositionalPhase phase : {
             MPMC::CompositionalPhase::Oil,
             MPMC::CompositionalPhase::Gas,
             MPMC::CompositionalPhase::Water})
    {
        MPMC::three_phase_detail::writeCompositionToPrimary<Indices>(
            primary, phaseComposition[static_cast<std::size_t>(MPMC::phaseIndex(phase))], phase);
    }

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence(static_cast<std::uint8_t>(6)); // G+W
    phaseState.phaseMoleFraction = {0.0, 0.70, 0.30};
    phaseState.overallComposition = {0.15, 0.45, 0.05, 0.35};

    const auto originalPrimary = primary;
    const auto originalState = phaseState;
    const auto result = equilibrium.updatePhaseState(primary, phaseState);

    require(result.recoverableFailure(),
            "invalid stability followed by failed unrestricted flash must be recoverable");
    require(result.stabilityInvalid,
            "failure diagnostics must record invalid stability explicitly");
    require(!result.missingPhaseUnstable,
            "invalid stability must not be mislabeled as certified missing-phase instability");
    require(result.unrestrictedFlashFailed,
            "invalid stability recovery must report unrestricted flash failure");
    require(primary == originalPrimary,
            "invalid-stability failure must restore primary state");
    require(phaseState.phasePresence.bits() == originalState.phasePresence.bits(),
            "invalid-stability failure must restore phase presence");
}


void testStableReducedPhaseSet()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    // Almost dry feed: ordinary PR settles on an oil+gas split with no
    // water-rich liquid.  This verifies that stability does not force all three
    // phases to exist merely because the model supports them.
    const Composition z{0.85, 0.05, 0.099999, 0.000001};
    const auto flash = equilibrium.flashPTZ(50.0e5, fluid.temperature, z);
    require(flash.converged, "reduced-phase PTz flash must converge");
    require(flash.presence.contains(MPMC::CompositionalPhase::Oil),
            "reduced state must contain oil-rich phase");
    require(flash.presence.contains(MPMC::CompositionalPhase::Gas),
            "reduced state must contain gas phase");
    require(!flash.presence.contains(MPMC::CompositionalPhase::Water),
            "stable dry state must keep the water-rich phase absent");
    near(flash.saturation[2], 0.0, 1.0e-14,
         "absent water-rich phase saturation must be zero");

    // Emulate a Newton iterate that still marks W active but drives Sw below
    // zero.  The post-removal O+G state is stable for this P-T-z, so the
    // restricted flash + stability state machine must keep W absent rather
    // than recreating it through an unrestricted flash.
    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 50.0e5;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, flash);
    phaseState.phasePresence.add(MPMC::CompositionalPhase::Water);
    phaseState.phaseMoleFraction[2] = 0.0;
    primary[Indices::Primary::waterSaturation] = -1.0e-6;

    equilibrium.updatePhaseState(primary, phaseState);
    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "stable reduced transition must retain oil-rich phase");
    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Gas),
            "stable reduced transition must retain gas phase");
    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Water),
            "stable missing water-rich phase must remain absent after reduction");
    near(primary[Indices::Primary::waterSaturation], 0.0, 1.0e-14,
         "stable absent water-rich phase saturation must remain zero");
}

} // namespace

int main()
{
    try
    {
        testThreePhasePTzFlash();
        testPhaseRemovalAndStabilityReappearance();
        testTraceAppearanceHoldPreventsImmediateDeletion();
        testFlowPhysicsUsesAllThreeCompositionalPhases();
        testRecoverablePhaseUpdateFailureIsTransactional();
        testInvalidStabilityFailureIsDistinguishedFromInstability();
        testStableReducedPhaseSet();
        std::cout << "Fully compositional PR three-phase flash: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Fully compositional PR three-phase flash: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
