/**
 * @file solve_result.hpp
 * @brief SNES/Newton/KSP 非线性求解结果与迭代历史数据结构。
 */
#pragma once

#include <string>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief PETSc/SNES 驱动记录的一次 Newton 迭代信息。
 *
 * `iteration == 0` is the initial nonlinear residual and therefore has no
 * associated linear solve.  For iteration N>0, `linearIterations` is the KSP
 * work used to advance from Newton state N-1 to N, and `wallTimeSeconds` is
 * the elapsed wall time for that Newton update.
 */
struct NewtonIterationRecord final
{
    int iteration{0};
    double residualNorm{0.0};
    double residualRatio{1.0};
    int linearIterations{0};
    long long cumulativeLinearIterations{0};
    double wallTimeSeconds{0.0};
    int linearReasonCode{0};
    std::string linearReason;
};

/** @brief 一次 SNES 求解摘要，包含 Newton/KSP 迭代历史。 */
struct NonlinearSolveResult final
{
    bool converged{false};
    int nonlinearIterations{0};
    int reasonCode{0};
    std::string reason;
    double wallTimeSeconds{0.0};

    long long linearIterations{0};
    double initialResidualNorm{0.0};
    double finalResidualNorm{0.0};
    double residualReduction{1.0};
    std::vector<NewtonIterationRecord> newtonHistory;
};

struct AdaptiveAttemptPlan final
{
    double timeStep{0.0};
    bool clippedToOutputTime{false};
};

} // namespace MPMC
