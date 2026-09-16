/**
 * @file callbacks.hpp
 * @brief PETSc SNES 残差/Jacobian 回调与 Natural 运行时桥接。
 */
#pragma once

#include <petscsnes.h>
#include <petscsys.h>

#include <exception>
#include <type_traits>
#include <utility>

namespace MPMC
{

namespace natural_callback_detail
{
template <class Runtime, class = void>
struct HasStateUpdateReport : std::false_type
{
};

template <class Runtime>
struct HasStateUpdateReport<
    Runtime,
    std::void_t<decltype(
        std::declval<Runtime &>().updateStateWithReport(std::declval<Vec>()))>>
    : std::true_type
{
};
} // namespace natural_callback_detail

/**
 * @brief 用 PETSc 标准 line-search pre-check 执行 Natural Newton 增量限制。
 *
 * 历史目标 PETSc 会在内部直接调用全局 `updateSol`。标准 PETSc 不包含该
 * 非公开扩展，因此必须把同一限制器注册到公开 pre-check 接口。重复调用是
 * 幂等的：已经满足限制的增量在第二次检查时保持不变。
 */
template <class Runtime>
PetscErrorCode naturalLineSearchPreCheck(
    SNESLineSearch,
    Vec solution,
    Vec step,
    PetscBool *changedStep,
    void *context)
{
    PetscFunctionBeginUser;

    const MPI_Comm comm = solution != nullptr
        ? PetscObjectComm(reinterpret_cast<PetscObject>(solution))
        : PETSC_COMM_SELF;
    PetscCheck(context != nullptr, comm, PETSC_ERR_ARG_NULL,
               "Natural line-search runtime context is null.");

    try
    {
        static_cast<Runtime *>(context)->limitNewtonStep(solution, step);
        // The limiter may have changed Y.  Conservatively report a change so
        // every PETSc line-search implementation rebuilds its trial state.
        if (changedStep != nullptr)
            *changedStep = PETSC_TRUE;
    }
    catch (const std::exception &error)
    {
        SETERRQ(comm, PETSC_ERR_LIB,
                "[ERROR][NEWTON] Natural line-search pre-check failed: %s",
                error.what());
    }
    catch (...)
    {
        SETERRQ(comm, PETSC_ERR_LIB,
                "[ERROR][NEWTON] Natural line-search pre-check failed with an unknown C++ exception.");
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief 用 PETSc 标准 line-search post-check 事务式更新 Natural 相态。
 *
 * phase-state 保存在独立 PETSc Vec 中，因此即使 trial solution 本身未被改写，
 * 也要把 `changedCandidate` 置真，保证 backtracking 等 line search 在相态更新后
 * 重新评价残差。
 */
template <class Runtime>
PetscErrorCode naturalLineSearchPostCheck(
    SNESLineSearch,
    Vec,
    Vec,
    Vec candidate,
    PetscBool *,
    PetscBool *changedCandidate,
    void *context)
{
    PetscFunctionBeginUser;

    const MPI_Comm comm = candidate != nullptr
        ? PetscObjectComm(reinterpret_cast<PetscObject>(candidate))
        : PETSC_COMM_SELF;
    PetscCheck(context != nullptr, comm, PETSC_ERR_ARG_NULL,
               "Natural line-search runtime context is null.");

    try
    {
        auto &runtime = *static_cast<Runtime *>(context);
        if constexpr (natural_callback_detail::HasStateUpdateReport<Runtime>::value)
        {
            const auto report = runtime.updateStateWithReport(candidate);
            if (!report.succeeded)
            {
                PetscCall(PetscInfo(
                    candidate,
                    "Natural post-check phase update failed: rank=%d cell=%lld "
                    "removed=%d unstable=%d stability_invalid=%d "
                    "restricted_fail=%d unrestricted_fail=%d\n",
                    report.firstFailureRank,
                    static_cast<long long>(report.firstFailureCell),
                    report.firstFailure.phaseRemoved ? 1 : 0,
                    report.firstFailure.missingPhaseUnstable ? 1 : 0,
                    report.firstFailure.stabilityInvalid ? 1 : 0,
                    report.firstFailure.restrictedFlashFailed ? 1 : 0,
                    report.firstFailure.unrestrictedFlashFailed ? 1 : 0));
            }
        }
        else
            runtime.updateState(candidate);

        if (changedCandidate != nullptr)
            *changedCandidate = PETSC_TRUE;
    }
    catch (const std::exception &error)
    {
        SETERRQ(comm, PETSC_ERR_LIB,
                "[ERROR][STATE ] Natural line-search post-check failed: %s",
                error.what());
    }
    catch (...)
    {
        SETERRQ(comm, PETSC_ERR_LIB,
                "[ERROR][STATE ] Natural line-search post-check failed with an unknown C++ exception.");
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief 将 Natural pre/post-check 安装到 SNES 当前实际使用的 line search。 */
template <class Runtime>
PetscErrorCode installNaturalLineSearchCallbacks(
    SNES snes,
    Runtime &runtime)
{
    PetscFunctionBeginUser;

    PetscCheck(snes != nullptr, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "SNES must not be null.");
    SNESLineSearch lineSearch = nullptr;
    PetscCall(SNESGetLineSearch(snes, &lineSearch));
    PetscCall(SNESLineSearchSetPreCheck(
        lineSearch,
        naturalLineSearchPreCheck<Runtime>,
        &runtime));
    PetscCall(SNESLineSearchSetPostCheck(
        lineSearch,
        naturalLineSearchPostCheck<Runtime>,
        &runtime));

    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief 标准 PETSc SNES residual 回调桥。
 *
 * Runtime 提供带 SNES 上下文的 residual 入口，从而可将局部热力学
 * 可恢复失败转换为 PETSc function-domain error，而不是 MPI_Abort。
 */
template <class Runtime>
PetscErrorCode naturalFormFunction(
    SNES snes,
    Vec solution,
    [[maybe_unused]] Vec residual,
    void *context)
{
    PetscFunctionBeginUser;

    [[maybe_unused]] const MPI_Comm comm =
        solution != nullptr
            ? PetscObjectComm(
                  reinterpret_cast<PetscObject>(
                      solution))
            : PETSC_COMM_SELF;

    PetscCheck(
        context != nullptr,
        comm,
        PETSC_ERR_ARG_NULL,
        "Natural SNES runtime context is null.");

    try
    {
        PetscCall(
            static_cast<Runtime *>(context)
                ->formFunction(
                    snes,
                    solution,
                    residual));
    }
    catch (const std::exception &error)
    {
        // 数值：Newton trial state 可能暂时离开 EOS 定义域。向 SNES 报告
        // function-domain failure，使外层时间步控制器可以拒绝并缩小 dt。
        PetscCall(VecSet(residual, 0.0));
        PetscCall(PetscInfo(snes, "Natural residual left its domain: %s\n", error.what()));
        PetscCall(SNESSetFunctionDomainError(snes));
    }
    catch (...)
    {
        PetscCall(VecSet(residual, 0.0));
        PetscCall(PetscInfo(snes, "Natural residual left its domain with an unknown C++ exception.\n"));
        PetscCall(SNESSetFunctionDomainError(snes));
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief 标准 PETSc SNES Jacobian 回调桥。
 */
template <class Runtime>
PetscErrorCode naturalFormJacobian(
    SNES,
    Vec solution,
    [[maybe_unused]] Mat jacobian,
    [[maybe_unused]] Mat preconditioner,
    void *context)
{
    PetscFunctionBeginUser;

    [[maybe_unused]] const MPI_Comm comm =
        solution != nullptr
            ? PetscObjectComm(
                  reinterpret_cast<PetscObject>(
                      solution))
            : PETSC_COMM_SELF;

    PetscCheck(
        context != nullptr,
        comm,
        PETSC_ERR_ARG_NULL,
        "Natural SNES runtime context is null.");

    try
    {
        PetscCall(
            static_cast<Runtime *>(context)
                ->formJacobian(
                    solution,
                    jacobian,
                    preconditioner));
    }
    catch (const std::exception &error)
    {
        SETERRQ(
            comm,
            PETSC_ERR_LIB,
            "[ERROR][JACOBIAN] Natural FormJacobian failed: %s",
            error.what());
    }
    catch (...)
    {
        SETERRQ(
            comm,
            PETSC_ERR_LIB,
            "[ERROR][JACOBIAN] Natural FormJacobian failed with an unknown C++ exception.");
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief 一次性把 Natural runtime 安装到标准 SNES callback。
 *
 * 标准 PETSc 通过 line-search pre/post-check 调用 Newton limiter 和相态更新；
 * 同时保留 SNES application context，兼容历史定制 PETSc 直接调用全局
 * `updateState/updateSol` 的路径。
 */
template <class Runtime>
PetscErrorCode installNaturalCallbacks(
    SNES snes,
    Runtime &runtime,
    Vec residual,
    Mat jacobian,
    Mat preconditioner = nullptr)
{
    PetscFunctionBeginUser;

    PetscCheck(
        snes != nullptr,
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_NULL,
        "SNES must not be null.");

    PetscCheck(
        residual != nullptr,
        PetscObjectComm(
            reinterpret_cast<PetscObject>(
                snes)),
        PETSC_ERR_ARG_NULL,
        "Residual Vec must not be null.");

    PetscCheck(
        jacobian != nullptr,
        PetscObjectComm(
            reinterpret_cast<PetscObject>(
                snes)),
        PETSC_ERR_ARG_NULL,
        "Jacobian Mat must not be null.");

    if (preconditioner == nullptr)
    {
        preconditioner = jacobian;
    }

    /*
     * 这里必须同时设置 SNES application context。
     *
     * 标准 SNESSetFunction()/SNESSetJacobian() 的 ctx 属于各自 callback；
     * 用户定制 PETSc 的 ls.c/virs.c/viss.c 则直接把 SNES application
     * context (`snes->user` in that customized source) 传给 updateState/updateSol。
     * 因此两类入口必须显式绑定同一个 Runtime。
     */
    PetscCall(
        SNESSetApplicationContext(
            snes,
            &runtime));

    PetscCall(
        SNESSetFunction(
            snes,
            residual,
            naturalFormFunction<Runtime>,
            &runtime));

    PetscCall(
        SNESSetJacobian(
            snes,
            jacobian,
            preconditioner,
            naturalFormJacobian<Runtime>,
            &runtime));

    PetscCall(installNaturalLineSearchCallbacks(snes, runtime));

    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief 用户定制 PETSc 调用的 `updateState(Vec,void*)` 的类型安全桥。
 *
 * 此函数本身不是全局符号。每个可执行程序只需定义一个很薄的全局包装：
 *
 * ```cpp
 * void updateState(Vec X, void *ctx)
 * {
 *     MPMC::naturalUpdateStateHook<MyRuntime>(X, ctx);
 * }
 * ```
 */
template <class Runtime>
void naturalUpdateStateHook(
    Vec solution,
    void *context) noexcept
{
    MPI_Comm comm = PETSC_COMM_WORLD;

    if (solution != nullptr)
    {
        comm = PetscObjectComm(
            reinterpret_cast<PetscObject>(
                solution));
    }

    if (context == nullptr)
    {
        PetscPrintf(
            comm,
            "[ERROR][STATE ] runtime context is null.\n");
        MPI_Abort(comm, PETSC_ERR_ARG_NULL);
        return;
    }

    try
    {
        auto &runtime = *static_cast<Runtime *>(context);
        if constexpr (natural_callback_detail::HasStateUpdateReport<Runtime>::value)
        {
            const auto report = runtime.updateStateWithReport(solution);
            if (!report.succeeded)
            {
                int rank = 0;
                MPI_Comm_rank(comm, &rank);
                if (rank == 0)
                {
                    PetscPrintf(
                        comm,
                        "[WARN ][STATE ] recoverable thermodynamic failure: "
                        "first_rank=%d cell=%lld failed_ranks=%zu removed=%d unstable=%d "
                        "stability_invalid=%d restricted_fail=%d unrestricted_fail=%d. "
                        "The state transaction was rolled back; the next SNES "
                        "residual evaluation will report a function-domain error.\n",
                        report.firstFailureRank,
                        static_cast<long long>(report.firstFailureCell),
                        report.failingRankCount,
                        report.firstFailure.phaseRemoved ? 1 : 0,
                        report.firstFailure.missingPhaseUnstable ? 1 : 0,
                        report.firstFailure.stabilityInvalid ? 1 : 0,
                        report.firstFailure.restrictedFlashFailed ? 1 : 0,
                        report.firstFailure.unrestrictedFlashFailed ? 1 : 0);
                }
            }
        }
        else
        {
            // Backward-compatible path for external/custom runtimes that only
            // implement the historical void updateState(Vec) contract.
            runtime.updateState(solution);
        }
    }
    catch (const std::exception &error)
    {
        PetscPrintf(
            comm,
            "[ERROR][STATE ] %s\n",
            error.what());
        MPI_Abort(comm, PETSC_ERR_LIB);
    }
    catch (...)
    {
        PetscPrintf(
            comm,
            "[ERROR][STATE ] unknown C++ exception.\n");
        MPI_Abort(comm, PETSC_ERR_LIB);
    }
}

/**
 * @brief 用户定制 PETSc 调用的 `updateSol(Vec,Vec,void*)` 的类型安全桥。
 *
 * 可执行程序应以与当前定制 PETSc 一致的 C++ linkage 定义全局 `updateSol` 包装。
 * Runtime::limitNewtonStep() 内部保持旧修改版 PETSc 的 Y 符号约定：
 * 进入时 Y 是线性求解得到的正 Newton correction，函数内部临时乘 -1、限制、
 * 再恢复符号，返回后 line-search 仍执行 `X <- X - lambda Y`。
 */
template <class Runtime>
void naturalUpdateSolHook(
    Vec solution,
    Vec step,
    void *context) noexcept
{
    MPI_Comm comm = PETSC_COMM_WORLD;

    if (solution != nullptr)
    {
        comm = PetscObjectComm(
            reinterpret_cast<PetscObject>(
                solution));
    }

    if (context == nullptr)
    {
        PetscPrintf(
            comm,
            "[ERROR][NEWTON] runtime context is null.\n");
        MPI_Abort(comm, PETSC_ERR_ARG_NULL);
        return;
    }

    try
    {
        static_cast<Runtime *>(context)
            ->limitNewtonStep(
                solution,
                step);
    }
    catch (const std::exception &error)
    {
        PetscPrintf(
            comm,
            "[ERROR][NEWTON] %s\n",
            error.what());
        MPI_Abort(comm, PETSC_ERR_LIB);
    }
    catch (...)
    {
        PetscPrintf(
            comm,
            "[ERROR][NEWTON] unknown C++ exception.\n");
        MPI_Abort(comm, PETSC_ERR_LIB);
    }
}

} // namespace MPMC
