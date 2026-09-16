/**
 * @file eos_parameter_regression_test.cpp
 * @brief 单元测试：验证 `eos_parameter_regression` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <tools/eos_parameter_regression.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using Config = MPMC::CompositionalModelConfig<
    2, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Composition = std::array<double, 2>;

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string &message)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (std::abs(actual - expected) > tolerance * scale)
        throw std::runtime_error(message);
}

Eos makeEos()
{
    constexpr std::array<double, 2> tc{305.32, 469.7};
    constexpr std::array<double, 2> pc{4.872e6, 3.37e6};
    constexpr std::array<double, 2> vc{1.45e-4, 3.13e-4};
    constexpr std::array<double, 2> omega{0.0995, 0.251};
    constexpr std::array<double, 2> mw{0.03007, 0.07215};
    constexpr std::array<std::array<double, 2>, 2> kij{{{{0.0, 0.0}}, {{0.0, 0.0}}}};
    return Eos(
        0.45724, 0.07780,
        MPMC::CompositionalMixture<Indices>(tc, pc, vc, omega, mw, kij),
        1, 2.414213562373095, -0.414213562373095, 1.0e-30);
}

void checkGenericOptimizer()
{
    const std::vector<MPMC::tools::RegressionParameter> parameters{
        {"a", 0.0, -5.0, 5.0, 0.5},
        {"b", 0.0, -5.0, 5.0, 0.5}};
    const auto result = MPMC::tools::boundedPatternSearch(
        parameters,
        [](const std::vector<double> &x) {
            return (x[0] - 1.25) * (x[0] - 1.25)
                + 2.0 * (x[1] + 0.75) * (x[1] + 0.75);
        });
    require(result.objective < 1.0e-9, "generic bounded pattern search must minimize a convex test objective");
    near(result.values[0], 1.25, 5.0e-5, "fitted quadratic parameter a");
    near(result.values[1], -0.75, 5.0e-5, "fitted quadratic parameter b");
}

void checkProductionEosBipRegression()
{
    constexpr double p = 4.0e6;
    constexpr double T = 350.0;
    const Composition x{0.45, 0.55};
    constexpr double targetKij = 0.1375;

    Eos target = makeEos();
    target.configureBinaryInteractionFunction(
        [](int i, int j, double) { return i == j ? 0.0 : targetKij; });
    const auto reference = target.phaseResult(
        p, T, x, MPMC::CompositionalPhase::Gas, false);

    const std::vector<MPMC::tools::RegressionParameter> parameters{
        {"k_C2_nC5", 0.0, -0.2, 0.4, 0.04}};
    MPMC::tools::PatternSearchOptions options;
    options.minimumRelativeStep = 1.0e-7;
    const auto result = MPMC::tools::boundedPatternSearch(
        parameters,
        [&](const std::vector<double> &candidate) {
            Eos eos = makeEos();
            const double kij = candidate[0];
            eos.configureBinaryInteractionFunction(
                [kij](int i, int j, double) { return i == j ? 0.0 : kij; });
            const auto phase = eos.phaseResult(
                p, T, x, MPMC::CompositionalPhase::Gas, false);
            double objective = 0.0;
            for (std::size_t c = 0; c < x.size(); ++c)
            {
                const double r = std::log(phase.fugacityCoefficient[c])
                    - std::log(reference.fugacityCoefficient[c]);
                objective += r * r;
            }
            return objective;
        },
        options);
    near(result.values[0], targetKij, 2.0e-5,
         "regression framework must recover a production PR BIP from EOS data");
    require(result.objective < result.initialObjective,
            "EOS regression objective must improve from its initial value");

    const auto path = std::filesystem::temp_directory_path() /
        "mpmc_eos_parameter_regression_test.csv";
    MPMC::tools::writeRegressionResultCsv(path, parameters, result);
    std::ifstream in(path);
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    require(text.find("k_C2_nC5") != std::string::npos,
            "regression CSV writer must persist parameter names");
    std::filesystem::remove(path);
}

} // namespace

int main()
{
    checkGenericOptimizer();
    checkProductionEosBipRegression();
    std::cout << "eos_parameter_regression_test: PASS\n";
    return 0;
}
