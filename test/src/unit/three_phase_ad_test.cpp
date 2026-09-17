/**
 * @file three_phase_ad_test.cpp
 * @brief 单元测试：验证 `three_phase_ad` 的核心语义、边界条件和回归行为。
 */
#include <common/math.hpp>
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

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
using Config = MPMC::CompositionalModelConfig<
    4, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ADIndices<Config>;
using Value = Indices::ValueType;
using Composition = std::array<double, Indices::numComponents>;

void require(bool condition, const std::string &message)
{
    if (!condition) throw std::runtime_error(message);
}

MPMC::FluidSystem<Indices> makeFluid()
{
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
    MPMC::CubicEquationOfState<Indices> eos(
        0.4572355, 0.0779691, std::move(mixture), 1,
        2.414213562373095, -0.414213562373095, 1.0e-30);
    MPMC::FluidSystem<Indices> fluid(
        {800.0, 20.0, 1000.0}, {1.0e-3, 1.0e-5, 1.0e-4}, std::move(eos),
        {"CO2", "CH4", "nC16", "H2O"}, {"Oil", "Gas", "Water"}, 350.0);
    fluid.configureFullyCompositionalThreePhase(3);
    fluid.gasRelativePermeability = [](Value s) { return s * s; };
    fluid.waterRelativePermeability = [](Value s) { return s * s; };
    fluid.threePhaseOilRelativePermeability = [](Value, Value so, Value) { return so * so; };
    return fluid;
}

void testAdPath()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    const Composition z{0.75, 0.025, 0.025, 0.20};
    const auto flash = equilibrium.flashPTZ(50.0e5, fluid.temperature, z);
    require(flash.converged && flash.presence.count() == 3,
            "AD test requires converged three-phase flash");

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = 50.0e5;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, flash);

    const auto state = MPMC::CellStateCodec<Indices>::decode(primary, phaseState);

    // Oil composition primaries are phase-component amounts q_i=S_o x_i.
    // The codec must recover the physical mole fractions while retaining both
    // the direct q derivative and the q/S_o chain rule.
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        const int qIndex = Indices::Primary::liquidComposition[c];
        const double so = flash.saturation[0];
        const double expectedQ = so * flash.composition[0][c];
        require(std::abs(primary[static_cast<std::size_t>(qIndex)] - expectedQ) < 1.0e-13,
                "PTz assignment must store Oil q_i=S_o x_i");
        require(std::abs(MPMC::scalarValue(state.liquidComponentAmount[c]) - expectedQ) < 1.0e-13,
                "CellState must retain the Oil phase-component amount");
        require(std::abs(MPMC::scalarValue(state.liquidMoleFraction[c]) - flash.composition[0][c]) < 1.0e-12,
                "Oil q/S_o decode must recover the physical mole fraction");
        require(std::abs(state.liquidComponentAmount[c].derivative(qIndex) - 1.0) < 1.0e-14,
                "Oil q coordinate must keep an identity AD column");
        require(std::abs(state.liquidMoleFraction[c].derivative(qIndex) - 1.0 / so) < 1.0e-9,
                "Oil mole fraction must expose d(q/S_o)/dq");
    }

    MPMC::CellPropertyEvaluator<Indices> evaluator(fluid);
    const auto properties = evaluator.evaluate(state, Value(0.25));
    const auto accumulation = MPMC::computeFluidAccumulation<Indices>(properties);

    static_assert(Value::size() == Indices::numPrimaryVariables,
                  "AD derivative dimension must follow the full primary layout");
    require(std::isfinite(MPMC::scalarValue(properties.density[0])),
            "AD oil density must be finite");
    require(std::isfinite(MPMC::scalarValue(properties.density[2])),
            "AD water-rich density must be finite");
    require(std::isfinite(MPMC::scalarValue(accumulation.componentMass[3])),
            "AD H2O accumulation must be finite");
    require(std::abs(properties.density[0].derivative(Indices::Primary::pressure)) > 0.0,
            "EOS density must retain pressure derivative in AD mode");

    auto exteriorPrimary = primary;
    exteriorPrimary[Indices::Primary::pressure] = 49.0e5;
    const auto exteriorState = MPMC::CellStateCodec<Indices>::decode(exteriorPrimary, phaseState);
    const auto exteriorProperties = evaluator.evaluate(exteriorState, Value(0.25));
    const auto flux = MPMC::computeFaceMassFlux<Indices>(
        state, properties, exteriorState, exteriorProperties,
        1.0e-12, 0.0, 1.0, 1.0, 0, 1);
    require(std::isfinite(MPMC::scalarValue(flux.component[3])),
            "AD three-phase H2O face flux must be finite");

    // G8K: residual uses a scalar-only face path.  Its values and recorded
    // upwind decisions must exactly match the AD kernel at the same state.
    MPMC::FaceFluxLinearizationDecision<Indices> faceDecision{};
    const auto valueFlux = MPMC::computeFaceMassFluxValue<Indices>(
        state, properties, exteriorState, exteriorProperties,
        1.0e-12, 0.0, 1.0, 1.0, 0, 1, -1, &faceDecision);
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        require(valueFlux.component[c] == MPMC::scalarValue(flux.component[c]),
                "G8K scalar face kernel changed a component flux value");
    }
    require(valueFlux.water == MPMC::scalarValue(flux.water),
            "G8K scalar face kernel changed water flux value");

    const auto decidedFlux = MPMC::computeFaceMassFlux<Indices>(
        state, properties, exteriorState, exteriorProperties,
        1.0e-12, 0.0, 1.0, 1.0, 0, 1, -1, &faceDecision);
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        require(MPMC::scalarValue(decidedFlux.component[c]) ==
                    MPMC::scalarValue(flux.component[c]),
                "G8K cached upwind decision changed face flux value");
        for (int derivative = 0; derivative < Indices::numPrimaryVariables; ++derivative)
        {
            require(decidedFlux.component[c].derivative(derivative) ==
                        flux.component[c].derivative(derivative),
                    "G8K cached upwind decision changed a face flux derivative");
        }
    }

    std::array<Value, Indices::numPhases> phasePressure{};
    std::array<Value, Indices::numPhases> density{};
    std::array<Value, Indices::numPhases> mobility{};
    for (int p = 0; p < Indices::numPhases; ++p)
    {
        const std::size_t i = static_cast<std::size_t>(p);
        phasePressure[i] = state.pressure;
        density[i] = properties.density[i];
        mobility[i] = properties.mobility[i];
    }
    std::array<double, Indices::numPhases> phaseFraction{};
    std::array<double, Indices::numComponents> injectionMass{};
    const auto well = MPMC::computePerforationWellSource<Indices>(
        MPMC::WellType::Producer, 1.0e-12, Value(45.0e5),
        phasePressure, density, mobility, fluid.surfaceDensity,
        phaseFraction, injectionMass, properties.massFraction);
    require(std::isfinite(MPMC::scalarValue(well.componentMassSource[3])),
            "AD three-phase H2O well source must be finite");

    // Instantiate the complete full-three-phase local residual as well.  This
    // catches equation-layout/template regressions that property/flux tests
    // alone cannot detect.
    const auto previousProperties = MPMC::scalarizeCellProperties<Indices>(properties);
    const auto residual = MPMC::assembleCellResidual<Indices>(
        state, properties, previousProperties, flux,
        well.componentMassSource, well.waterMassSource,
        1.0, 86400.0, 1.0);
    for (const auto &equation : residual.value)
        require(std::isfinite(MPMC::scalarValue(equation)),
                "AD full-three-phase local residual must be finite");

    // G8J: accepted previous accumulation may be cached once and reused by
    // residual/Jacobian.  The split local/base + transport/source path must
    // remain exactly equivalent in both values and AD derivatives.
    const auto previousAccumulation =
        MPMC::computeFluidAccumulation<Indices>(previousProperties);
    auto splitResidual = MPMC::assembleCellLocalResidual<Indices>(
        state, properties, previousProperties,
        1.0, 86400.0, 1.0, -1, 2650.0, nullptr,
        &previousAccumulation, nullptr);
    MPMC::addCellConservationTransportAndWell<Indices>(
        splitResidual, flux, well.componentMassSource, well.waterMassSource);

    for (int equation = 0; equation < Indices::numEquations; ++equation)
    {
        const auto &expected = residual.value[static_cast<std::size_t>(equation)];
        const auto &actual = splitResidual.value[static_cast<std::size_t>(equation)];
        require(MPMC::scalarValue(actual) == MPMC::scalarValue(expected),
                "G8J split residual changed an equation value");
        for (int derivative = 0; derivative < Indices::numPrimaryVariables; ++derivative)
        {
            require(actual.derivative(derivative) == expected.derivative(derivative),
                    "G8J split residual changed an AD derivative");
        }
    }

}

void testOilMinorityContinuation()
{
    auto fluid = makeFluid();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    const Composition z{0.75, 0.025, 0.025, 0.20};
    const auto flash = equilibrium.flashPTZ(50.0e5, fluid.temperature, z);
    require(flash.converged && flash.presence.count() == 3,
            "minority continuation regression requires a three-phase flash");

    std::array<double, Indices::numPrimaryVariables> base{};
    base[Indices::Primary::pressure] = 50.0e5;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(base, phaseState, flash);

    // Perturb Gas composition so the O-G fugacity block is intentionally off
    // equilibrium, then hold the physical O/G/W compositions fixed while only
    // changing S_o. The continuation row must scale linearly with S_o.
    base[static_cast<std::size_t>(Indices::Primary::vaporComposition[0])] += 1.0e-4;

    const auto makeMinorityState = [&](double so) {
        auto primary = base;
        primary[Indices::Primary::liquidSaturation] = so;
        primary[Indices::Primary::vaporSaturation] =
            1.0 - so - primary[Indices::Primary::waterSaturation];
        for (int component = 0;
             component < Indices::numIndependentCompositionsPerPhase;
             ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            primary[static_cast<std::size_t>(Indices::Primary::liquidComposition[c])] =
                so * flash.composition[0][c];
        }
        return primary;
    };

    MPMC::CellPropertyEvaluator<Indices> evaluator(fluid);
    MPMC::CellProperties<Indices, double> previous{};
    const auto residualAt = [&](double so) {
        const auto primary = makeMinorityState(so);
        const auto state = MPMC::CellStateCodec<Indices>::decode(primary, phaseState);
        const auto properties = evaluator.evaluate(state, Value(0.25));
        return MPMC::assembleCellLocalResidual<Indices>(
            state, properties, previous, 1.0, 1.0, 1.0);
    };

    constexpr double largeSo = 1.0e-2;
    constexpr double smallSo = 1.0e-3;
    const auto largeResidual = residualAt(largeSo);
    const auto smallResidual = residualAt(smallSo);
    const auto equation = static_cast<std::size_t>(Indices::Equation::fugacity[0]);
    const double largeValue = MPMC::scalarValue(largeResidual.value[equation]);
    const double smallValue = MPMC::scalarValue(smallResidual.value[equation]);
    require(std::abs(smallValue) > 1.0e-16,
            "continuation regression needs a nonzero O-G mismatch");
    require(std::abs(largeValue / smallValue - largeSo / smallSo) < 2.0e-8,
            "O-* fugacity residual must vanish linearly with S_o at fixed composition");

    // At the inactive endpoint, the same equation block changes to q_i=0 and
    // q_N=S_o-sum(q_i)=0.  Check the exact AD closure, including dependent q_N.
    auto inactivePrimary = makeMinorityState(0.0);
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        inactivePrimary[static_cast<std::size_t>(
            Indices::Primary::liquidComposition[static_cast<std::size_t>(component)])] = 0.0;
    }
    auto inactivePhaseState = phaseState;
    inactivePhaseState.phasePresence.remove(MPMC::CompositionalPhase::Oil);
    const auto inactiveState =
        MPMC::CellStateCodec<Indices>::decode(inactivePrimary, inactivePhaseState);
    MPMC::CellProperties<Indices, Value> inactiveProperties{};
    const auto inactiveResidual = MPMC::assembleCellLocalResidual<Indices>(
        inactiveState, inactiveProperties, previous, 1.0, 1.0, 1.0);

    const int q0 = Indices::Primary::liquidComposition[0];
    require(std::abs(inactiveResidual.value[static_cast<std::size_t>(
                         Indices::Equation::fugacity[0])].derivative(q0) - 1.0) < 1.0e-14,
            "inactive Oil q_i row must keep an identity Jacobian");
    const auto &dependentRow = inactiveResidual.value[
        static_cast<std::size_t>(Indices::Equation::fugacity.back())];
    require(std::abs(dependentRow.derivative(Indices::Primary::liquidSaturation) - 1.0) < 1.0e-14,
            "inactive Oil q_N row must contain +dS_o");
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        const int qIndex = Indices::Primary::liquidComposition[static_cast<std::size_t>(component)];
        require(std::abs(dependentRow.derivative(qIndex) + 1.0) < 1.0e-14,
                "inactive Oil q_N row must contain -dq_i");
    }
}

} // namespace

int main()
{
    try
    {
        testAdPath();
        testOilMinorityContinuation();
        std::cout << "Fully compositional three-phase AD path: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Fully compositional three-phase AD path: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
