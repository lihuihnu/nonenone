/**
 * @file cpa_single_nonaqueous_phase_test.cpp
 * @brief 回归二维 CPA 算例中富 nC10 单非水相被误标为 Gas 的问题。
 */
#include "../../../case/h2o_co2_nc10_2d_benchmark/case_config.hpp"
#include "../../../case/h2o_co2_nc10_2d_benchmark/case_fluid.hpp"

#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/state/three_phase_equilibrium.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    false,
    false,
    false,
    false,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;
using Composition = std::array<double, Indices::numComponents>;
using PrimaryArray = std::array<double, Indices::numPrimaryVariables>;
using Equilibrium = MPMC::FullyCompositionalThreePhaseEquilibrium<Indices>;

struct BenchmarkState final
{
    int inputIndex;
    double pressure;
    double gasSaturation;
    double waterSaturation;
    Composition gasComposition;
    Composition waterComposition;
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

Composition composition(double x0, double x1)
{
    return Composition{x0, x1, 1.0 - x0 - x1};
}

void writeComposition(
    PrimaryArray& primary,
    const std::array<int, Indices::numIndependentCompositionsPerPhase>& indices,
    const Composition& x)
{
    for (int component = 0;
         component < Indices::numIndependentCompositionsPerPhase;
         ++component)
    {
        primary[static_cast<std::size_t>(
            indices[static_cast<std::size_t>(component)])] =
            x[static_cast<std::size_t>(component)];
    }
}

void checkCpaSingleNonaqueousState(const BenchmarkState& state)
{
    auto fluid =
        MPMC::cases::makeFluidSystem<Indices, CaseConfig::CpaFactoryConfig>();
    BenchmarkCaseFluid::apply(fluid);
    const Equilibrium equilibrium(fluid);

    PrimaryArray primary{};
    primary[Indices::Primary::pressure] = state.pressure;
    primary[Indices::Primary::liquidSaturation] = 0.0;
    primary[Indices::Primary::vaporSaturation] = state.gasSaturation;
    primary[Indices::Primary::waterSaturation] = state.waterSaturation;
    writeComposition(
        primary,
        Indices::Primary::liquidComposition,
        state.gasComposition);
    writeComposition(
        primary,
        Indices::Primary::vaporComposition,
        state.gasComposition);
    writeComposition(
        primary,
        Indices::Primary::waterComposition,
        state.waterComposition);

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence::gasOnly();
    phaseState.phasePresence.add(MPMC::CompositionalPhase::Water);
    equilibrium.updateSecondary(primary, phaseState);

    const auto overallBefore = phaseState.overallComposition;
    equilibrium.updatePhaseState(primary, phaseState);

    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "CPA dense single nonaqueous phase must be relabeled Oil for input_index " +
                std::to_string(state.inputIndex));
    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Gas),
            "CPA dense single nonaqueous phase must not remain Gas for input_index " +
                std::to_string(state.inputIndex));
    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Water),
            "CPA reduced state must retain Water for input_index " +
                std::to_string(state.inputIndex));
    require(
        primary[Indices::Primary::liquidSaturation] > 0.75,
        "CPA relabeled oil saturation must carry the former dense nonaqueous phase");
    require(primary[Indices::Primary::vaporSaturation] == 0.0,
            "CPA relabeled state must zero inactive gas saturation");
    for (std::size_t component = 0; component < overallBefore.size(); ++component)
    {
        require(
            std::abs(phaseState.overallComposition[component] -
                     overallBefore[component]) < 2.0e-12,
            "CPA single nonaqueous relabeling must preserve overall composition");
    }
}
} // namespace

int main()
{
    try {
        const std::array<BenchmarkState, 3> states{{
            {308,
             4.82065490572822746e6,
             7.96985725006652435e-1,
             2.03014274993347482e-1,
             composition(2.63180684194782178e-3,
                         2.07197603712508310e-1),
             composition(9.93258583537875439e-1,
                         6.74141625275343197e-3)},
            {905,
             4.82395152294778544e6,
             7.96034448270393935e-1,
             2.03965551729606120e-1,
             composition(2.86661845445388948e-3,
                         2.77909685400717088e-1),
             composition(9.91026416340109684e-1,
                         8.97358344300016111e-3)},
            {121,
             4.82724732571661193e6,
             7.94805800087118097e-1,
             2.05194199912881819e-1,
             composition(3.11968245897279199e-3,
                         3.42997679787385212e-1),
             composition(9.89046176773337171e-1,
                         1.09538230044786723e-2)}
        }};

        for (const auto& state : states)
            checkCpaSingleNonaqueousState(state);

        std::cout
            << "CPA single nonaqueous phase canonicalization: ALL PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "CPA single nonaqueous phase canonicalization: FAIL: "
            << error.what() << '\n';
        return 1;
    }
}
