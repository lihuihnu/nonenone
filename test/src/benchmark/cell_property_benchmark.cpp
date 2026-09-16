/**
 * @file cell_property_benchmark.cpp
 * @brief 性能基准：测量 `cell_property_benchmark` 对应核心路径的计算开销。
 */
#include "../../../case/3p4c_pr_reservoir/case_config.hpp"

#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/state/cell_property_evaluator.hpp>
#include <natural/state/state_codec.hpp>
#include <natural/state/three_phase_equilibrium.hpp>

#include <array>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>

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
} // namespace

int main()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::Config>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    const auto flash = equilibrium.flashPTZ(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);
    if (!flash.converged)
        return 2;

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = CaseConfig::InitialState::pressure;
    MPMC::PhaseStateData<Indices> phaseState{};
    equilibrium.assignFlashResult(primary, phaseState, flash);

    const auto state = MPMC::CellStateCodec<Indices>::decode(primary, phaseState);
    MPMC::CellPropertyEvaluator<Indices> evaluator(fluid);

    constexpr std::size_t iterations = 50000;
    double checksum = 0.0;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < iterations; ++i)
    {
        const auto properties = evaluator.evaluate(state, CaseConfig::Rock::porosity);
        checksum += properties.density[0] +
            properties.viscosity[1] +
            properties.fugacity[2][0];
    }
    const auto stop = std::chrono::steady_clock::now();

    const double elapsedMs =
        std::chrono::duration<double, std::milli>(stop - start).count();
    const double microsecondsPerEvaluation =
        elapsedMs * 1000.0 / static_cast<double>(iterations);

    std::cout << std::fixed << std::setprecision(3)
              << "iterations=" << iterations
              << " elapsed_ms=" << elapsedMs
              << " us_per_eval=" << microsecondsPerEvaluation
              << " checksum=" << std::setprecision(15) << checksum
              << '\n';
    return 0;
}
