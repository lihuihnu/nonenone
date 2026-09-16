/**
 * @file eos_parameter_regression.hpp
 * @brief EOS 经验参数的离线有界回归工具。
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace MPMC::tools
{

/** @brief 一个在 production solver 外部拟合的有界经验参数。 */
struct RegressionParameter
{
    std::string name{};
    double initialValue{0.0};
    double lowerBound{0.0};
    double upperBound{0.0};
    double initialStep{0.0};
};

struct PatternSearchOptions
{
    int maximumIterations{80};
    double contractionFactor{0.5};
    double expansionFactor{1.2};
    double minimumRelativeStep{1.0e-6};
    double improvementTolerance{1.0e-12};
};

struct RegressionResult
{
    bool converged{false};
    int iterations{0};
    int objectiveEvaluations{0};
    double initialObjective{std::numeric_limits<double>::quiet_NaN()};
    double objective{std::numeric_limits<double>::quiet_NaN()};
    std::vector<double> values{};
    std::vector<double> finalSteps{};
};

namespace regression_detail
{

inline void validateParameter(const RegressionParameter &parameter)
{
    if (parameter.name.empty())
        throw std::invalid_argument("Regression parameter name cannot be empty.");
    if (!std::isfinite(parameter.initialValue) ||
        !std::isfinite(parameter.lowerBound) ||
        !std::isfinite(parameter.upperBound) ||
        !std::isfinite(parameter.initialStep))
    {
        throw std::invalid_argument("Regression parameter values must be finite.");
    }
    if (!(parameter.lowerBound < parameter.upperBound))
        throw std::invalid_argument("Regression parameter lower bound must be below upper bound.");
    if (parameter.initialValue < parameter.lowerBound ||
        parameter.initialValue > parameter.upperBound)
    {
        throw std::invalid_argument("Regression parameter initial value must lie inside its bounds.");
    }
    if (parameter.initialStep < 0.0)
        throw std::invalid_argument("Regression parameter initial step cannot be negative.");
}

inline double defaultStep(const RegressionParameter &parameter)
{
    return parameter.initialStep > 0.0
        ? parameter.initialStep
        : 0.1 * (parameter.upperBound - parameter.lowerBound);
}

inline double relativeStep(
    double step,
    const RegressionParameter &parameter)
{
    const double scale = std::max({
        1.0,
        std::abs(parameter.lowerBound),
        std::abs(parameter.upperBound),
        std::abs(parameter.initialValue)});
    return step / scale;
}

} // namespace regression_detail

/**
 * @brief 确定性的无导数有界坐标/模式搜索算法。
 *
 * The optimizer deliberately lives in tools/.  Production PR/SW/CPA code only
 * exposes physical parameters; an objective callback reconstructs/configures
 * the selected EOS and evaluates experimental or synthetic data.  This keeps
 * parameter fitting completely outside the reservoir nonlinear solve.
 */
template <class Objective>
[[nodiscard]] RegressionResult boundedPatternSearch(
    const std::vector<RegressionParameter> &parameters,
    Objective &&objective,
    PatternSearchOptions options = {})
{
    if (parameters.empty())
        throw std::invalid_argument("EOS regression requires at least one parameter.");
    if (options.maximumIterations <= 0)
        throw std::invalid_argument("Regression maximum iteration count must be positive.");
    if (!(options.contractionFactor > 0.0 && options.contractionFactor < 1.0) ||
        !(options.expansionFactor >= 1.0) ||
        !(options.minimumRelativeStep > 0.0) ||
        !(options.improvementTolerance >= 0.0) ||
        !std::isfinite(options.contractionFactor) ||
        !std::isfinite(options.expansionFactor) ||
        !std::isfinite(options.minimumRelativeStep) ||
        !std::isfinite(options.improvementTolerance))
    {
        throw std::invalid_argument("Invalid pattern-search options.");
    }

    std::vector<double> x;
    std::vector<double> step;
    x.reserve(parameters.size());
    step.reserve(parameters.size());
    for (const auto &parameter : parameters)
    {
        regression_detail::validateParameter(parameter);
        x.push_back(parameter.initialValue);
        step.push_back(regression_detail::defaultStep(parameter));
    }

    RegressionResult result;
    auto evaluate = [&](const std::vector<double> &candidate) {
        const double value = objective(candidate);
        ++result.objectiveEvaluations;
        if (!std::isfinite(value))
            return std::numeric_limits<double>::infinity();
        return value;
    };

    double best = evaluate(x);
    if (!std::isfinite(best))
        throw std::runtime_error("Initial EOS regression objective is not finite.");
    result.initialObjective = best;

    for (int iteration = 0; iteration < options.maximumIterations; ++iteration)
    {
        bool improvedAny = false;
        for (std::size_t i = 0; i < parameters.size(); ++i)
        {
            const auto &parameter = parameters[i];
            std::vector<double> bestCandidate = x;
            double localBest = best;

            for (double direction : {-1.0, 1.0})
            {
                auto candidate = x;
                candidate[i] = std::clamp(
                    x[i] + direction * step[i],
                    parameter.lowerBound,
                    parameter.upperBound);
                if (candidate[i] == x[i])
                    continue;
                const double candidateObjective = evaluate(candidate);
                const double threshold = options.improvementTolerance *
                    std::max(1.0, std::abs(localBest));
                if (candidateObjective + threshold < localBest)
                {
                    localBest = candidateObjective;
                    bestCandidate = std::move(candidate);
                }
            }

            if (localBest < best)
            {
                x = std::move(bestCandidate);
                best = localBest;
                step[i] *= options.expansionFactor;
                step[i] = std::min(
                    step[i], parameter.upperBound - parameter.lowerBound);
                improvedAny = true;
            }
        }

        result.iterations = iteration + 1;
        if (!improvedAny)
        {
            bool allSmall = true;
            for (std::size_t i = 0; i < step.size(); ++i)
            {
                step[i] *= options.contractionFactor;
                if (regression_detail::relativeStep(step[i], parameters[i]) >
                    options.minimumRelativeStep)
                {
                    allSmall = false;
                }
            }
            if (allSmall)
            {
                result.converged = true;
                break;
            }
        }
    }

    result.values = x;
    result.finalSteps = step;
    result.objective = best;
    if (!result.converged)
    {
        result.converged = true;
        for (std::size_t i = 0; i < step.size(); ++i)
        {
            if (regression_detail::relativeStep(step[i], parameters[i]) >
                options.minimumRelativeStep)
            {
                result.converged = false;
                break;
            }
        }
    }
    return result;
}

/** @brief 在不绑定特定 EOS 参数模式的前提下保存拟合值。 */
inline void writeRegressionResultCsv(
    const std::filesystem::path &path,
    const std::vector<RegressionParameter> &parameters,
    const RegressionResult &result)
{
    if (parameters.size() != result.values.size() ||
        parameters.size() != result.finalSteps.size())
    {
        throw std::invalid_argument("Regression result size does not match parameter descriptors.");
    }
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create EOS regression CSV: " + path.string());
    out << std::setprecision(17)
        << "parameter,initial,lower,upper,fitted,final_step\n";
    for (std::size_t i = 0; i < parameters.size(); ++i)
    {
        out << parameters[i].name << ','
            << parameters[i].initialValue << ','
            << parameters[i].lowerBound << ','
            << parameters[i].upperBound << ','
            << result.values[i] << ','
            << result.finalSteps[i] << '\n';
    }
    out << "objective_initial," << result.initialObjective << ",,,objective_final,"
        << result.objective << '\n';
    out << "converged," << (result.converged ? 1 : 0)
        << ",,,iterations," << result.iterations << '\n';
    out << "objective_evaluations," << result.objectiveEvaluations << ",,,,\n";
}

} // namespace MPMC::tools
