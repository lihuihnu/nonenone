/**
 * @file natural_conservation_test.cpp
 * @brief 单元测试：验证 `natural_conservation` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <natural/fluid_system.hpp>
#include <natural/physics/face_flux.hpp>
#include <natural/state/cell_property_evaluator.hpp>
#include <natural/state/cell_state.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

using Config =
    MPMC::CompositionalModelConfig<
        6, true, true,
        false, false, false>;
using Indices = MPMC::ScalarIndices<Config>;

void near(double a, double b, double tol, const std::string &message)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    if (std::abs(a - b) > tol * scale)
        throw std::runtime_error(message);
}

MPMC::CompositionalMixture<Indices> mixture()
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

MPMC::FluidSystem<Indices> fluid()
{
    MPMC::CubicEquationOfState<Indices> eos(
        0.4572355,
        0.0779691,
        mixture(),
        1,
        2.414213562373095,
        -0.414213562373095);

    MPMC::FluidSystem<Indices> fs(
        {800.0, 2.0, 1000.0},
        {1.0e-4, 2.0e-5, 2.0e-4},
        std::move(eos),
        {"N2/CH4", "CO2", "C2-5", "C6-13", "C14-24", "C25-80"},
        {"liquid", "vapor", "water"},
        387.45);

    fs.gasRelativePermeability = [](double s) { return s * s; };
    fs.oilRelativePermeability = [](double s) { return s * s; };
    fs.waterRelativePermeability = [](double s) { return s * s; };
    fs.waterViscosity = [](double) { return 2.0e-4; };
    fs.waterFormationVolumeFactor = [](double) { return 1.0; };
    fs.threePhaseOilRelativePermeability = [](double, double so, double) { return so * so; };
    return fs;
}

MPMC::CellState<Indices, double> state(double pressure)
{
    MPMC::CellState<Indices, double> s{};
    s.pressure = pressure;
    s.liquidSaturation = 0.55;
    s.vaporSaturation = 0.25;
    s.waterSaturation = 0.20;
    s.liquidMoleFraction = {
        0.3246914, 0.0128351, 0.2278401, 0.2606985, 0.1134144, 0.0605205};
    s.vaporMoleFraction = {
        0.808671, 0.025284, 0.1487798, 0.0175878, 0.0006759, 0.0000015};
    s.hydrocarbonPhaseState = MPMC::HydrocarbonPhaseState::TwoPhase;
    return s;
}

} // namespace

int main()
{
    auto fs = fluid();
    MPMC::CellPropertyEvaluator<Indices> evaluator(fs);

    const auto leftState = state(15.5e6);
    const auto rightState = state(14.5e6);
    const auto left = evaluator.evaluate(leftState, 0.25);
    const auto right = evaluator.evaluate(rightState, 0.25);

    constexpr double transmissibility = 2.0e-12;
    constexpr double leftVolume = 1000.0;
    constexpr double rightVolume = 900.0;

    const auto leftOut = MPMC::computeFaceMassFlux<Indices>(
        leftState, left, rightState, right,
        transmissibility, 0.0,
        leftVolume, rightVolume,
        0, 1);

    const auto rightOut = MPMC::computeFaceMassFlux<Indices>(
        rightState, right, leftState, left,
        transmissibility, 0.0,
        rightVolume, leftVolume,
        1, 0);

    for (int component = 0; component < Indices::numComponents; ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        near(
            leftOut.component[c] + rightOut.component[c],
            0.0,
            2.0e-12,
            "internal-face component mass must cancel pairwise");
    }

    near(
        leftOut.water + rightOut.water,
        0.0,
        2.0e-12,
        "internal-face water mass must cancel pairwise");

    std::cout << "Natural two-cell conservation: ALL PASS\n";
    return 0;
}
