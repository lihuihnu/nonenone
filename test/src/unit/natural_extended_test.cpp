/**
 * @file natural_extended_test.cpp
 * @brief 单元测试：验证 `natural_extended` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <natural/assembly/cell_residual.hpp>
#include <natural/fluid_system.hpp>
#include <natural/kernel/cell_kernel.hpp>
#include <natural/physics/face_flux.hpp>
#include <natural/state/cell_property_evaluator.hpp>
#include <natural/state/phase_equilibrium.hpp>
#include <natural/state/state_codec.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace
{

using BaseConfig =
    MPMC::CompositionalModelConfig<
        6, true, true,
        false, false, false>;
using Base = MPMC::ScalarIndices<BaseConfig>;
using AdBase = MPMC::ADIndices<BaseConfig>;

using FullConfig =
    MPMC::CompositionalModelConfig<
        6, true, true,
        true, true, true>;
using AdFull = MPMC::ADIndices<FullConfig>;

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(
    double actual,
    double expected,
    double tolerance,
    const std::string &message)
{
    const double scale =
        std::max({1.0, std::abs(actual), std::abs(expected)});

    if (std::abs(actual - expected) >
        tolerance * scale)
    {
        throw std::runtime_error(
            message +
            ": actual=" + std::to_string(actual) +
            ", expected=" + std::to_string(expected));
    }
}

template <class Indices>
MPMC::CompositionalMixture<Indices> makeMixture()
{
    return MPMC::CompositionalMixture<Indices>(
        {189.515, 304.2, 387.607, 597.497, 698.515, 875.0},
        {4580011.59, 7386592.50, 4095515.97,
         3345244.875, 1768374.5625, 1169006.79},
        {9.97012032965401e-5, 9.26344713533338e-5,
         0.000217076707259486, 0.000381162235869935,
         0.000721410148917871, 0.00113570073874421},
        {0.00854, 0.228, 0.16733, 0.38609, 0.80784, 1.23141},
        {0.0161594, 0.04401, 0.0455725, 0.11774, 0.248827, 0.48152},
        {{0.0, 0.00070981, 0.00077754, 0.0100, 0.0110, 0.0110},
         {0.00070981, 0.0, 0.1500, 0.1500, 0.1500, 0.1500},
         {0.00077754, 0.1500, 0.0, 0.0, 0.0, 0.0},
         {0.0100, 0.1500, 0.0, 0.0, 0.0, 0.0},
         {0.0110, 0.1500, 0.0, 0.0, 0.0, 0.0},
         {0.0110, 0.1500, 0.0, 0.0, 0.0, 0.0}});
}

template <class Indices>
auto makeEos()
{
    return MPMC::CubicEquationOfState<Indices>(
        0.4572355,
        0.0779691,
        makeMixture<Indices>(),
        1,
        2.414213562373095,
        -0.414213562373095);
}

template <class Indices>
auto makeFluid()
{
    using Fluid = MPMC::FluidSystem<Indices>;
    using Value = typename Indices::ValueType;

    Fluid fluid(
        {800.0, 2.0, 1000.0},
        {1.0e-4, 2.0e-5, 2.0e-4},
        makeEos<Indices>(),
        {"N2/CH4", "CO2", "C2-5", "C6-13", "C14-24", "C25-80"},
        {"liquid", "vapor", "water"},
        387.45);

    fluid.gasRelativePermeability = [](Value s)
    {
        return s * s;
    };

    fluid.oilRelativePermeability = [](Value s)
    {
        return s * s;
    };

    fluid.waterRelativePermeability = [](Value s)
    {
        return s * s;
    };

    fluid.waterViscosity = [](Value)
    {
        return Value(2.0e-4);
    };

    fluid.waterFormationVolumeFactor = [](Value)
    {
        return Value(1.0);
    };

    fluid.threePhaseOilRelativePermeability = [](
        Value,
        Value so,
        Value)
    {
        return so * so;
    };

    return fluid;
}

std::array<double, AdFull::numPrimaryVariables>
makeFullPrimary()
{
    std::array<double, AdFull::numPrimaryVariables> primary{};

    primary[AdFull::Primary::pressure] = 15.0e6;

    const std::array<double, 5> x{
        0.3246914,
        0.0128351,
        0.2278401,
        0.2606985,
        0.1134144};

    const std::array<double, 5> y{
        0.808671,
        0.025284,
        0.1487798,
        0.0175878,
        0.0006759};

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        primary[AdFull::Primary::liquidComposition[i]] = x[i];
        primary[AdFull::Primary::vaporComposition[i]] = y[i];
    }

    primary[AdFull::Primary::waterSaturation] = 0.2;
    primary[AdFull::Primary::wellPressure] = 10.0e6;
    primary[AdFull::Primary::liquidSaturation] = 0.5;
    primary[AdFull::Primary::vaporSaturation] = 0.3;
    primary[AdFull::Primary::aqueousCO2MoleFraction] = 0.01;

    return primary;
}

void testCellCodecAndPropertyEvaluator()
{
    auto fluid = makeFluid<AdFull>();

    fluid.configureAqueousCO2(
        1,
        0.01801528,
        0.0);

    fluid.configureAdsorption(
        {1.468e-3, 2.860e-3, 0.0, 0.0, 0.0, 0.0},
        {1.471e-7, 2.472e-7, 0.0, 0.0, 0.0, 0.0});

    const auto primary = makeFullPrimary();

    MPMC::PhaseStateData<AdFull> phaseState;
    phaseState.phase = MPMC::HydrocarbonPhaseState::TwoPhase;
    phaseState.overallComposition = {
        0.4630, 0.0164, 0.2052, 0.19108, 0.08113, 0.04319};

    const auto state =
        MPMC::CellStateCodec<AdFull>::decode(
            primary,
            phaseState);

    near(
        state.liquidMoleFraction.back().value(),
        0.0605205,
        1.0e-12,
        "dependent liquid composition");

    near(
        state.vaporMoleFraction.back().value(),
        -0.0009985,
        1.0e-12,
        "dependent vapor composition");

    MPMC::CellPropertyEvaluator<AdFull> evaluator(fluid);

    std::array<AdFull::ValueType, AdFull::numComponents> z{};
    for (int c = 0; c < AdFull::numComponents; ++c)
        z[static_cast<std::size_t>(c)] =
            phaseState.overallComposition[static_cast<std::size_t>(c)];

    const auto properties = evaluator.evaluate(
        state,
        AdFull::ValueType(0.25),
        0.6,
        2.0,
        &z);

    for (int phase = 0; phase < AdFull::numPhases; ++phase)
    {
        require(
            std::isfinite(
                properties.density[
                    static_cast<std::size_t>(phase)].value()),
            "cell density must be finite");

        require(
            std::isfinite(
                properties.mobility[
                    static_cast<std::size_t>(phase)].value()),
            "cell mobility must be finite");
    }

    require(
        properties.trappedGasSaturation.value() > 0.0,
        "Land branch should trap gas during imbibition");

    require(
        properties.aqueousCO2MassFraction.value() > 0.0,
        "aqueous CO2 mass fraction should be positive");

    require(
        properties.adsorbedVolume[0].value() > 0.0 &&
        properties.adsorbedVolume[1].value() > 0.0,
        "competitive adsorption should produce positive loading");

    require(
        std::isfinite(
            properties.density[AdFull::Phase::liquid].derivative(
                AdFull::Primary::pressure)),
        "cell property AD pressure derivative must be finite");

    const double liquidMassSum =
        std::accumulate(
            properties.massFraction[AdFull::Phase::liquid].begin(),
            properties.massFraction[AdFull::Phase::liquid].end(),
            0.0,
            [](double sum, const auto &value)
            {
                return sum + value.value();
            });

    near(
        liquidMassSum,
        1.0,
        1.0e-12,
        "liquid mass fractions close");
}

void testPhaseStateUpdate()
{
    auto fluid = makeFluid<Base>();

    std::array<double, Base::numPrimaryVariables> primary{};
    primary[Base::Primary::pressure] = 15.0e6;

    const std::array<double, 5> x{
        0.3246914,
        0.0128351,
        0.2278401,
        0.2606985,
        0.1134144};

    const std::array<double, 5> y{
        0.808671,
        0.025284,
        0.1487798,
        0.0175878,
        0.0006759};

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        primary[Base::Primary::liquidComposition[i]] = x[i];
        primary[Base::Primary::vaporComposition[i]] = y[i];
    }

    primary[Base::Primary::waterSaturation] = 0.2;
    primary[Base::Primary::wellPressure] = 10.0e6;
    primary[Base::Primary::liquidSaturation] = 0.5;
    primary[Base::Primary::vaporSaturation] = 0.3;

    MPMC::PhaseStateData<Base> phaseState;
    phaseState.phase = MPMC::HydrocarbonPhaseState::TwoPhase;

    MPMC::PhaseEquilibriumManager<Base> manager(fluid);
    manager.updateSecondary(primary, phaseState);

    const double zSum =
        std::accumulate(
            phaseState.overallComposition.begin(),
            phaseState.overallComposition.end(),
            0.0);

    near(zSum, 1.0, 1.0e-12, "phase-state overall composition closure");
    require(
        std::isfinite(phaseState.liquidCompressibility) &&
        std::isfinite(phaseState.vaporCompressibility),
        "phase-state Z factors must be finite");

    // 状态预处理只做物理上下界裁剪，不额外把三相饱和度归一化写回。
    auto boundedPrimary = primary;
    boundedPrimary[Base::Primary::liquidSaturation] = 0.8;
    boundedPrimary[Base::Primary::vaporSaturation] = 0.4;
    boundedPrimary[Base::Primary::waterSaturation] = 0.2;
    manager.sanitizePrimaryBeforeFlash(
        boundedPrimary,
        false);

    near(
        boundedPrimary[Base::Primary::liquidSaturation],
        0.8,
        0.0,
        "phase sanitize must not rewrite saturation normalization");
}

void testExteriorLinearization()
{
    using Eval = AdBase::ValueType;

    MPMC::CellState<AdBase, Eval> interiorState{};
    MPMC::CellState<AdBase, Eval> exteriorState{};
    MPMC::CellProperties<AdBase, Eval> interior{};
    MPMC::CellProperties<AdBase, Eval> exterior{};

    interiorState.pressure =
        Eval::createVariable(
            100.0,
            AdBase::Primary::pressure);

    exteriorState.pressure =
        Eval::createVariable(
            90.0,
            AdBase::Primary::pressure);

    interior.density[AdBase::Phase::liquid] = 100.0;
    exterior.density[AdBase::Phase::liquid] = 100.0;
    interior.saturation[AdBase::Phase::liquid] = 1.0;
    exterior.saturation[AdBase::Phase::liquid] = 1.0;

    const auto difference =
        MPMC::phasePotentialDifference<AdBase, Eval>(
            interiorState,
            interior,
            exteriorState,
            exterior,
            AdBase::Phase::liquid,
            0.0,
            1.0,
            1.0,
            0,
            1);

    near(
        difference.first.derivative(
            AdBase::Primary::pressure),
        -1.0,
        0.0,
        "exterior pressure must be AD-frozen in a local cell assembly");
}

void testVaporOnlyDegenerateResidual()
{
    using Eval = AdBase::ValueType;

    std::array<double, AdBase::numPrimaryVariables> primary{};
    primary[AdBase::Primary::pressure] = 1.0e7;
    primary[AdBase::Primary::waterSaturation] = 0.2;
    primary[AdBase::Primary::wellPressure] = 9.0e6;
    primary[AdBase::Primary::liquidSaturation] = 0.0;
    primary[AdBase::Primary::vaporSaturation] = 0.8;

    for (std::size_t c = 0;
         c < AdBase::Primary::liquidComposition.size();
         ++c)
    {
        primary[AdBase::Primary::liquidComposition[c]] =
            c == 0 ? 0.2 : 0.1;
        primary[AdBase::Primary::vaporComposition[c]] =
            c == 0 ? 0.3 : 0.1;
    }

    auto state =
        MPMC::CellStateCodec<AdBase>::decode(
            primary,
            MPMC::HydrocarbonPhaseState::VaporOnly);

    MPMC::CellProperties<AdBase, Eval> properties{};
    MPMC::CellProperties<AdBase, double> oldProperties{};
    MPMC::FaceMassFlux<AdBase, Eval> flux{};
    std::array<Eval, AdBase::numComponents> well{};

    const auto residual =
        MPMC::assembleCellResidual<AdBase, Eval>(
            state,
            properties,
            oldProperties,
            flux,
            well,
            Eval(0.0),
            1.0,
            1.0,
            1.0);

    const int row =
        AdBase::Equation::fugacity[0];

    near(
        residual.value[row].derivative(
            AdBase::Primary::vaporComposition[0]),
        1.0,
        0.0,
        "vapor-only degenerate row must constrain vapor composition");

    near(
        residual.value[row].derivative(
            AdBase::Primary::liquidComposition[0]),
        0.0,
        0.0,
        "vapor-only degenerate row must not constrain liquid composition");
}


void testAdIndicesDoubleKFunctionDecay()
{
    auto fluid = makeFluid<AdBase>();
    using Eval = typename AdBase::ValueType;

    typename decltype(fluid.eos)::EquilibriumConstantFunctions kFunctions{};
    for (int component = 0; component < AdBase::numComponents; ++component)
    {
        const double offset = 0.01 * static_cast<double>(component);
        kFunctions[static_cast<std::size_t>(component)] =
            [offset](Eval, Eval, const std::array<Eval, AdBase::numComponents> &z)
            {
                return Eval(1.0 + offset) + 0.5 * z[0];
            };
    }
    fluid.eos.configureEquilibriumConstants(std::move(kFunctions));

    const std::array<double, AdBase::numComponents> z{
        0.30, 0.10, 0.20, 0.15, 0.15, 0.10};

    const auto K =
        fluid.eos.evaluateEquilibriumConstants(
            15.0e6,
            387.45,
            z);

    static_assert(
        std::is_same_v<
            typename decltype(K)::value_type,
            double>,
        "double phase-state flash must receive scalar K values");

    near(
        K[0],
        1.15,
        1.0e-14,
        "ADIndices K-function must decay to double in scalar flash");

    near(
        K[5],
        1.20,
        1.0e-14,
        "ADIndices K-function scalar decay must preserve component value");
}

void testBlackOilAdIndicesWithScalarState()
{
    // Production PETSc runtimes use ADIndices, but initialization/diagnostic
    // paths may evaluate secondary properties from plain double states.  The
    // black-oil callbacks therefore return AD values internally while this
    // public call must decay them safely back to double.
    auto fluid = makeFluid<AdBase>();

    using Eval = typename AdBase::ValueType;
    fluid.blackOilProperties.configure(
        {[](const Eval &, const Eval &) { return Eval(700.0); },
         [](const Eval &, const Eval &) { return Eval(100.0); }},
        {[](const Eval &, const Eval &) { return Eval(1.0e-4); },
         [](const Eval &, const Eval &) { return Eval(2.0e-5); }});

    const std::array<double, AdBase::numComponents> composition{
        0.2, 0.2, 0.2, 0.15, 0.15, 0.1};

    near(
        fluid.blackOilProperties.computeMolarDensity(
            1.0e7, composition, true),
        700.0,
        1.0e-14,
        "ADIndices black-oil density must decay safely to scalar double");

    near(
        fluid.blackOilProperties.computeViscosity(
            1.0e7, composition, false),
        2.0e-5,
        1.0e-14,
        "ADIndices black-oil viscosity must decay safely to scalar double");
}

void testBlackOilOverallCompositionPath()
{
    auto fluid = makeFluid<Base>();

    typename decltype(fluid.eos)::EquilibriumConstantFunctions kFunctions{};
    for (int component = 0; component < Base::numComponents; ++component)
    {
        const double offset = 0.01 * component;
        kFunctions[static_cast<std::size_t>(component)] =
            [offset](double, double, const std::array<double, Base::numComponents> &z)
            {
                return 1.0 + z[0] + offset;
            };
    }
    fluid.eos.configureEquilibriumConstants(std::move(kFunctions));

    fluid.blackOilProperties.configure(
        {[](const double &, const double &) { return 700.0; },
         [](const double &, const double &) { return 100.0; }},
        {[](const double &, const double &) { return 1.0e-4; },
         [](const double &, const double &) { return 2.0e-5; }});

    MPMC::CellState<Base, double> state{};
    state.pressure = 1.0e7;
    state.liquidSaturation = 0.6;
    state.vaporSaturation = 0.2;
    state.waterSaturation = 0.2;
    state.hydrocarbonPhaseState = MPMC::HydrocarbonPhaseState::TwoPhase;
    state.liquidMoleFraction = {0.2, 0.2, 0.2, 0.15, 0.15, 0.1};
    state.vaporMoleFraction = {0.4, 0.2, 0.15, 0.1, 0.1, 0.05};

    const std::array<double, Base::numComponents> overall{
        0.7, 0.1, 0.05, 0.05, 0.05, 0.05};

    MPMC::CellPropertyEvaluator<Base> evaluator(fluid);
    const auto properties = evaluator.evaluate(
        state,
        0.25,
        0.0,
        2.0,
        &overall);

    const double expectedK0 = 1.0 + overall[0];
    near(
        properties.fugacity[0][0],
        expectedK0 * state.liquidMoleFraction[0],
        1.0e-14,
        "black-oil K path must use phase-state overall composition");
}

void testDeterministicPhaseDisappearance()
{
    auto fluid = makeFluid<Base>();
    MPMC::PhaseEquilibriumManager<Base> manager(fluid);

    std::array<double, Base::numPrimaryVariables> primary{};
    primary[Base::Primary::pressure] = 15.0e6;
    primary[Base::Primary::waterSaturation] = 0.2;
    primary[Base::Primary::wellPressure] = 10.0e6;
    primary[Base::Primary::liquidSaturation] = 0.8;
    primary[Base::Primary::vaporSaturation] = 0.0;

    for (std::size_t c = 0; c < Base::Primary::liquidComposition.size(); ++c)
    {
        primary[Base::Primary::liquidComposition[c]] = c == 0 ? 0.4 : 0.1;
        primary[Base::Primary::vaporComposition[c]] = c == 0 ? 0.4 : 0.1;
    }

    MPMC::PhaseStateData<Base> phaseState;
    phaseState.phase = MPMC::HydrocarbonPhaseState::TwoPhase;
    phaseState.liquidMoleFraction = 1.0;
    phaseState.equilibriumRatio.fill(1.0);

    manager.updatePhaseState(primary, phaseState);

    require(
        phaseState.phase == MPMC::HydrocarbonPhaseState::LiquidOnly,
        "zero vapor saturation must switch a two-phase cell to liquid-only");
    near(
        primary[Base::Primary::liquidSaturation],
        0.8,
        1.0e-14,
        "liquid-only saturation preserves non-water pore fraction");
    near(
        primary[Base::Primary::vaporSaturation],
        0.0,
        0.0,
        "vapor saturation disappears");
}


void testSinglePhasePreservesEquilibriumRatio()
{
    auto fluid = makeFluid<Base>();
    MPMC::PhaseEquilibriumManager<Base> manager(fluid);

    std::array<double, Base::numPrimaryVariables> primary{};
    primary[Base::Primary::pressure] = 15.0e6;
    primary[Base::Primary::waterSaturation] = 0.2;
    primary[Base::Primary::liquidSaturation] = 0.8;
    primary[Base::Primary::vaporSaturation] = 0.0;

    const std::array<double, Base::numComponents> composition{
        0.463, 0.0164, 0.2052, 0.19108, 0.08113, 0.04319};
    for (int c = 0; c < Base::numIndependentCompositionsPerPhase; ++c)
    {
        const auto i = static_cast<std::size_t>(c);
        primary[Base::Primary::liquidComposition[i]] = composition[i];
        // In a single-phase state this is only an inactive finite placeholder.
        primary[Base::Primary::vaporComposition[i]] = composition[i];
    }

    const std::array<double, Base::numComponents> originalK{
        2.42902000269171,
        1.93302307385877,
        0.646223364847938,
        0.0683987270072828,
        0.00663493246768742,
        5.15345236540681e-05};

    MPMC::PhaseStateData<Base> phaseState;
    phaseState.phase = MPMC::HydrocarbonPhaseState::LiquidOnly;
    phaseState.liquidMoleFraction = 1.0;
    phaseState.equilibriumRatio = originalK;
    phaseState.overallComposition = composition;

    manager.updateSecondary(primary, phaseState);

    for (int c = 0; c < Base::numComponents; ++c)
    {
        const auto i = static_cast<std::size_t>(c);
        near(
            phaseState.equilibriumRatio[i],
            originalK[i],
            1.0e-14,
            "single-phase secondary update must preserve historical K");
    }
}

void testLiquidOnlyCanCreateVaporPhase()
{
    auto fluid = makeFluid<Base>();
    MPMC::PhaseEquilibriumManager<Base> manager(fluid);

    // This 3p6c feed is unstable as a single liquid phase for the same PR EOS
    // at 150 bar.  It is therefore a compact regression for phase appearance.
    const std::array<double, Base::numComponents> z{
        0.463, 0.0164, 0.2052, 0.19108, 0.08113, 0.04319};
    const std::array<double, Base::numComponents> K{
        2.42902000269171,
        1.93302307385877,
        0.646223364847938,
        0.0683987270072828,
        0.00663493246768742,
        5.15345236540681e-05};

    const auto stability = manager.testStability(z, K, 15.0e6);
    require(!stability.stable, "phase-appearance regression feed must be unstable");

    std::array<double, Base::numPrimaryVariables> primary{};
    primary[Base::Primary::pressure] = 15.0e6;
    primary[Base::Primary::waterSaturation] = 0.2;
    primary[Base::Primary::liquidSaturation] = 0.8;
    primary[Base::Primary::vaporSaturation] = 0.0;
    for (int c = 0; c < Base::numIndependentCompositionsPerPhase; ++c)
    {
        const auto i = static_cast<std::size_t>(c);
        primary[Base::Primary::liquidComposition[i]] = z[i];
        primary[Base::Primary::vaporComposition[i]] = z[i];
    }

    MPMC::PhaseStateData<Base> phaseState;
    phaseState.phase = MPMC::HydrocarbonPhaseState::LiquidOnly;
    phaseState.liquidMoleFraction = 1.0;
    phaseState.equilibriumRatio = K;
    phaseState.overallComposition = z;

    manager.updatePhaseState(primary, phaseState);

    require(
        phaseState.phase == MPMC::HydrocarbonPhaseState::TwoPhase,
        "unstable liquid-only cell must create the vapor phase");
    require(
        primary[Base::Primary::vaporSaturation] > 0.0,
        "newly created vapor phase must receive positive seed saturation");
    require(
        primary[Base::Primary::liquidSaturation] > 0.0,
        "phase appearance must retain a positive liquid saturation");
    for (int c = 0; c < Base::numIndependentCompositionsPerPhase; ++c)
    {
        const auto i = static_cast<std::size_t>(c);
        near(primary[Base::Primary::liquidComposition[i]], z[i], 1.0e-14,
             "vapor appearance must preserve the existing liquid composition");
    }
}

void testNaturalCellKernelFacade()
{
    auto fluid = makeFluid<AdFull>();
    fluid.configureAqueousCO2(1, 0.01801528);
    fluid.configureAdsorption(
        {1.468e-3, 2.860e-3, 0, 0, 0, 0},
        {1.471e-7, 2.472e-7, 0, 0, 0, 0});

    MPMC::NaturalCellKernel<AdFull> kernel(fluid);

    const auto primary = makeFullPrimary();
    MPMC::PhaseStateData<AdFull> phaseState;
    phaseState.phase = MPMC::HydrocarbonPhaseState::TwoPhase;
    phaseState.overallComposition = {
        0.4630, 0.0164, 0.2052, 0.19108, 0.08113, 0.04319};

    const auto state = kernel.decodeState(primary, phaseState);

    std::array<AdFull::ValueType, AdFull::numComponents> z{};
    for (int c = 0; c < AdFull::numComponents; ++c)
        z[static_cast<std::size_t>(c)] =
            phaseState.overallComposition[static_cast<std::size_t>(c)];

    const auto properties =
        kernel.evaluateProperties(
            state,
            AdFull::ValueType(0.25),
            0.6,
            2.0,
            &z);

    require(
        std::isfinite(
            properties.porosity.value()),
        "NaturalCellKernel facade property result finite");
}

} // namespace

int main()
{
    testCellCodecAndPropertyEvaluator();
    testPhaseStateUpdate();
    testExteriorLinearization();
    testVaporOnlyDegenerateResidual();
    testAdIndicesDoubleKFunctionDecay();
    testBlackOilAdIndicesWithScalarState();
    testBlackOilOverallCompositionPath();
    testDeterministicPhaseDisappearance();
    testSinglePhasePreservesEquilibriumRatio();
    testLiquidOnlyCanCreateVaporPhase();
    testNaturalCellKernelFacade();

    std::cout
        << "Natural extended regression: ALL PASS\n";

    return 0;
}
