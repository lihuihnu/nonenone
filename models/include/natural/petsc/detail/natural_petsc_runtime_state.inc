/**
 * @file natural_petsc_runtime_state.inc
 * @brief Natural/PETSc 相态初始化、事务更新、历史提交和 Newton 限制实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    /**
     * @brief 用用户给定的统一 secondary state 初始化全部本地单元。
     *
     * 该入口用于严格复现已有基准算例：旧算例往往直接给出相态 flag、
     * K 值、总体组成 z 和液相摩尔分率 L，而不是从当前 EOS 重新 flash。
     * 函数只写 phase-state，不修改 primary solution；调用后仍应执行
     * initializeHistory(solution) 保存首个时间层。
     */
    void initializeUniformPhaseState(const PhaseStateData<Indices> &initialState)
    {
        const auto &phaseMap =
            backend_.dofMap(static_cast<PetscInt>(Indices::numPhaseStateVariables));

        PetscOwnedWriteView phaseView(phaseState_);
        for (const Cell cell : ownedCells_)
        {
            writeOwnedPhaseState_(
                phaseView.data(),
                phaseView.begin(),
                phaseView.end(),
                phaseMap,
                cell,
                initialState);
        }
        phaseStateInitialized_ = true;
    }

    /**
     * @brief 用当前解初始化 secondary phase-state，但不执行相切换。
     *
     * 推荐流程：
     * `fill X -> initializePhaseStateFromSolution -> updateState -> initializeHistory`。
     */
    void initializePhaseStateFromSolution(
        Vec solution,
        HydrocarbonPhaseState initialPhase =
            HydrocarbonPhaseState::TwoPhase)
    {
        static_assert(!Indices::fullyCompositionalThreePhase,
                      "FullyCompositionalThreePhase must be initialized from P-T-z with initializeUniformFromPTZ().");
        requireVector_(solution, "solution");

        const auto &primaryMap =
            backend_.dofMap(
                static_cast<PetscInt>(
                    Indices::numPrimaryVariables));

        const auto &phaseMap =
            backend_.dofMap(
                static_cast<PetscInt>(
                    Indices::numPhaseStateVariables));

        PetscOwnedReadView solutionView(
            solution);

        PetscOwnedWriteView phaseView(
            phaseState_);

        for (const Cell cell :
             ownedCells_)
        {
            PrimaryArray primary =
                readOwnedPrimary_(
                    solutionView.data(),
                    solutionView.begin(),
                    solutionView.end(),
                    primaryMap,
                    cell);

            PhaseStateData<Indices> state;
            state.phase = initialPhase;

            kernel_.phaseEquilibrium()
                .updateSecondary(
                    primary,
                    state);

            writeOwnedPhaseState_(
                phaseView.data(),
                phaseView.begin(),
                phaseView.end(),
                phaseMap,
                cell,
                state);
        }

        phaseStateInitialized_ = true;
    }

    /**
     * @brief 仅由 P-T-z 初始化全组分模型的所有单元。
     *
     * The current Natural flow model is isothermal, therefore @p temperature
     * must match FluidSystem::temperature.  The PR three-phase flash determines
     * the active phase set, phase mole fractions, all three phase compositions,
     * compressibility factors and pore-volume saturations.  The caller only
     * supplies pressure [Pa], temperature [K] and normalized/normalizable overall
     * mole fractions.
     */
    void initializeUniformFromPTZ(
        Vec solution,
        double pressure,
        double temperature,
        const std::array<double, Indices::numComponents> &overallComposition)
    {
        static_assert(Indices::fullyCompositionalThreePhase,
                      "initializeUniformFromPTZ is only available for FullyCompositionalThreePhase.");
        requireVector_(solution, "solution");
        if (!(pressure > 0.0) || !std::isfinite(pressure))
            throw std::invalid_argument("PTz initialization pressure must be finite and positive.");
        if (!(temperature > 0.0) || !std::isfinite(temperature))
            throw std::invalid_argument("PTz initialization temperature must be finite and positive.");
        const double temperatureTolerance =
            1.0e-10 * std::max(1.0, std::abs(fluid_.temperature));
        if (std::abs(temperature - fluid_.temperature) > temperatureTolerance)
            throw std::invalid_argument(
                "Natural runtime is isothermal: PTz initialization temperature must match FluidSystem::temperature.");

        const auto flash = kernel_.phaseEquilibrium().flashPTZ(
            pressure, temperature, overallComposition);
        if (!flash.converged)
            throw std::runtime_error("Uniform PTz three-phase flash did not converge.");

        const auto &primaryMap = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables));
        const auto &phaseMap = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPhaseStateVariables));

        PetscOwnedWriteView solutionView(solution);
        PetscOwnedWriteView phaseView(phaseState_);
        for (const Cell cell : ownedCells_)
        {
            PrimaryArray primary{};
            primary[Indices::Primary::pressure] = pressure;
            if constexpr (Indices::hasWellUnknown)
                primary[Indices::Primary::wellPressure] = pressure;

            PhaseStateData<Indices> state;
            kernel_.phaseEquilibrium().assignFlashResult(primary, state, flash);

            writeOwnedPrimary_(
                solutionView.data(), solutionView.begin(), solutionView.end(),
                primaryMap, cell, primary);
            writeOwnedPhaseState_(
                phaseView.data(), phaseView.begin(), phaseView.end(),
                phaseMap, cell, state);
        }
        phaseStateInitialized_ = true;
    }

    /**
     * @brief 事务式更新相态；热力学失败时不提交任何 cell 的部分更新。
     *
     * fully-compositional 相态管理器可以返回可恢复热力学失败。Runtime 先在
     * C++ scratch 中完成所有 owned cell 的状态转换，只有全部成功后才一次性
     * 写回 PETSc Vec。因此某个 cell 的 restricted/unrestricted flash 失败不会
     * 留下“前半部分 cell 已更新、后半部分仍是旧状态”的混合快照。
     */
    StateUpdateReport updateStateWithReport(Vec solution)
    {
        requirePhaseStateInitialized_();
        requireVector_(solution, "solution");

        PetscObjectState solutionStateBefore{};
        PetscObjectState phaseStateBefore{};
        PetscCallAbort(
            backend_.communicator(),
            VecGetState(solution, &solutionStateBefore));
        PetscCallAbort(
            backend_.communicator(),
            VecGetState(phaseState_, &phaseStateBefore));

        if (phaseEquilibriumCacheValid_ &&
            phaseEquilibriumSolution_ == solution &&
            phaseEquilibriumSolutionState_ == solutionStateBefore &&
            phaseEquilibriumPhaseStateState_ == phaseStateBefore)
        {
            return lastStateUpdateReport_;
        }

        const auto &primaryMap = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables));
        const auto &phaseMap = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPhaseStateVariables));

        stateUpdateScratch_.clear();
        if (stateUpdateScratch_.capacity() < ownedCells_.size())
            stateUpdateScratch_.reserve(ownedCells_.size());

        StateUpdateReport report{};

        {
            PetscOwnedReadView solutionView(solution);
            PetscOwnedReadView phaseView(phaseState_);

            for (const Cell cell : ownedCells_)
            {
                PendingStateUpdate pending;
                pending.cell = cell;
                pending.primary = readOwnedPrimary_(
                    solutionView.data(), solutionView.begin(), solutionView.end(),
                    primaryMap, cell);
                pending.phaseState = readOwnedPhaseState_(
                    phaseView.data(), phaseView.begin(), phaseView.end(),
                    phaseMap, cell);

                kernel_.phaseEquilibrium().sanitizePrimaryBeforeFlash(
                    pending.primary, options_.useVariableBounds);

                const PhaseUpdateResult update =
                    kernel_.phaseEquilibrium().updatePhaseState(
                        pending.primary, pending.phaseState);

                if (update.recoverableFailure())
                {
                    report.succeeded = false;
                    report.failingRankCount = 1;
                    report.firstFailureCell = cell;
                    report.firstFailure = update;
                    break;
                }

                ++report.localUpdatedCellCount;
                if (update.status == PhaseUpdateStatus::StableReducedSet)
                    ++report.localStableReducedCellCount;
                stateUpdateScratch_.push_back(std::move(pending));
            }
        }

        // updateState is an MPI transaction: if one rank cannot construct a
        // certified local thermodynamic state, no rank commits its scratch.
        PetscInt localFailure = report.succeeded ? 0 : 1;
        PetscInt globalFailure = 0;
        PetscCallAbort(
            backend_.communicator(),
            MPI_Allreduce(
                &localFailure, &globalFailure, 1, MPIU_INT, MPI_MAX,
                backend_.communicator()));

        if (globalFailure != 0)
        {
            int rank = 0;
            PetscCallAbort(
                backend_.communicator(),
                MPI_Comm_rank(backend_.communicator(), &rank));

            const int localFailureRank = report.succeeded
                ? std::numeric_limits<int>::max()
                : rank;
            int firstFailureRank = localFailureRank;
            PetscCallAbort(
                backend_.communicator(),
                MPI_Allreduce(
                    &localFailureRank, &firstFailureRank, 1, MPI_INT, MPI_MIN,
                    backend_.communicator()));

            PetscInt localFailureCount =
                static_cast<PetscInt>(report.failingRankCount);
            PetscInt globalFailureCount = 0;
            PetscCallAbort(
                backend_.communicator(),
                MPI_Allreduce(
                    &localFailureCount, &globalFailureCount, 1, MPIU_INT,
                    MPI_SUM, backend_.communicator()));

            PetscInt failureCell =
                rank == firstFailureRank
                    ? static_cast<PetscInt>(report.firstFailureCell)
                    : PetscInt(-1);
            std::array<PetscInt, 6> failurePayload{};
            if (rank == firstFailureRank)
            {
                failurePayload = {
                    static_cast<PetscInt>(report.firstFailure.status),
                    report.firstFailure.phaseRemoved ? PetscInt(1) : PetscInt(0),
                    report.firstFailure.missingPhaseUnstable ? PetscInt(1) : PetscInt(0),
                    report.firstFailure.stabilityInvalid ? PetscInt(1) : PetscInt(0),
                    report.firstFailure.restrictedFlashFailed ? PetscInt(1) : PetscInt(0),
                    report.firstFailure.unrestrictedFlashFailed ? PetscInt(1) : PetscInt(0)};
            }
            PetscCallAbort(
                backend_.communicator(),
                MPI_Bcast(
                    &failureCell, 1, MPIU_INT, firstFailureRank,
                    backend_.communicator()));
            PetscCallAbort(
                backend_.communicator(),
                MPI_Bcast(
                    failurePayload.data(),
                    static_cast<int>(failurePayload.size()), MPIU_INT,
                    firstFailureRank, backend_.communicator()));

            report.succeeded = false;
            report.failingRankCount =
                static_cast<std::size_t>(globalFailureCount);
            report.firstFailureRank = firstFailureRank;
            report.firstFailureCell = static_cast<Cell>(failureCell);
            report.firstFailure.status =
                static_cast<PhaseUpdateStatus>(failurePayload[0]);
            report.firstFailure.phaseRemoved = failurePayload[1] != 0;
            report.firstFailure.missingPhaseUnstable = failurePayload[2] != 0;
            report.firstFailure.stabilityInvalid = failurePayload[3] != 0;
            report.firstFailure.restrictedFlashFailed = failurePayload[4] != 0;
            report.firstFailure.unrestrictedFlashFailed = failurePayload[5] != 0;

            stateUpdateScratch_.clear();
            lastStateUpdateReport_ = report;
            pendingThermoStateFailure_ = true;
            phaseEquilibriumCacheValid_ = false;
            currentCacheValid_ = false;
            return report;
        }

        {
            PetscOwnedWriteView solutionView(solution);
            PetscOwnedWriteView phaseView(phaseState_);
            for (const auto &pending : stateUpdateScratch_)
            {
                writeOwnedPrimary_(
                    solutionView.data(), solutionView.begin(), solutionView.end(),
                    primaryMap, pending.cell, pending.primary);
                writeOwnedPhaseState_(
                    phaseView.data(), phaseView.begin(), phaseView.end(),
                    phaseMap, pending.cell, pending.phaseState);
            }
        }
        stateUpdateScratch_.clear();

        PetscObjectState solutionStateAfter{};
        PetscObjectState phaseStateAfter{};
        PetscCallAbort(
            backend_.communicator(),
            VecGetState(solution, &solutionStateAfter));
        PetscCallAbort(
            backend_.communicator(),
            VecGetState(phaseState_, &phaseStateAfter));

        phaseEquilibriumSolution_ = solution;
        phaseEquilibriumSolutionState_ = solutionStateAfter;
        phaseEquilibriumPhaseStateState_ = phaseStateAfter;
        phaseEquilibriumCacheValid_ = true;
        currentCacheValid_ = false;
        pendingThermoStateFailure_ = false;
        lastStateUpdateReport_ = report;
        return report;
    }

    /**
     * @brief 兼容既有 case/custom-PETSc hook 的状态更新入口。
     *
     * 保持 v67.8 及更早版本的 `void updateState(Vec)` 调用契约；需要结构化
     * 失败信息的新代码应调用 `updateStateWithReport()` 或读取
     * `lastStateUpdateReport()`。
     */
    void updateState(Vec solution)
    {
        (void)updateStateWithReport(solution);
    }

    /**
     * @brief 回滚失败时间步的解、相态以及仅属于该 Newton 尝试的瞬态状态。
     *
     * A recoverable phase-update failure sets pendingThermoStateFailure_.  A
     * timestep retry must clear that latch after restoring the last accepted
     * vectors; otherwise every smaller retry is rejected before evaluating its
     * own state.
     */
    void rollbackTimeStepAttempt(Vec solution)
    {
        requirePhaseStateInitialized_();
        requireVector_(solution, "solution");

        PetscCallAbort(
            backend_.communicator(),
            VecCopy(previousSolution_, solution));
        PetscCallAbort(
            backend_.communicator(),
            VecCopy(previousPhaseState_, phaseState_));

        stateUpdateScratch_.clear();
        lastStateUpdateReport_ = {};
        pendingThermoStateFailure_ = false;
        phaseEquilibriumCacheValid_ = false;
        currentCacheValid_ = false;
    }

    /** @brief 最近一次 updateState 的结构化状态；用于 timestep/诊断层决定是否重试。 */
    [[nodiscard]] const StateUpdateReport &lastStateUpdateReport() const noexcept
    {
        return lastStateUpdateReport_;
    }

    /**
     * @brief 保存一个已收敛时间层，供下一时间步蓄积项使用。
     *
     * 同时更新 Land 历史最大气饱和度。不要在每个 Newton 迭代调用。
     */
    void commitTimeStep(Vec solution)
    {
        requirePhaseStateInitialized_();
        requireVector_(solution, "solution");

        const bool reuseCurrentStateCache =
            canReuseCurrentCacheForHistory_(solution);

        PetscCallAbort(
            backend_.communicator(),
            VecCopy(
                solution,
                previousSolution_));

        PetscCallAbort(
            backend_.communicator(),
            VecCopy(
                phaseState_,
                previousPhaseState_));

        updateMaximumGasSaturation_(
            solution);

        // Land history only changes on accepted steps. Scatter it once here and
        // let all Newton residual/Jacobian evaluations reuse the local snapshot.
        refreshMaximumGasLocalCache_();

        if (reuseCurrentStateCache)
            refreshPreviousPropertiesFromCurrentCache_();
        else
            refreshPreviousPropertiesCache_();

        currentCacheValid_ = false;
        historyInitialized_ = true;
    }

    /** 首个时间层初始化与 commitTimeStep 同义。 */
    void initializeHistory(Vec solution)
    {
        commitTimeStep(solution);
    }

    /**
     * @brief 限制 SNES Newton step，保持原 updateSol 的符号约定。
     *
     * 原定制 PETSc 传入的 Y 先取负成为 `deltaSol`，限制后再取负写回。
     */
    void limitNewtonStep(
        Vec solution,
        Vec step)
    {
        requirePhaseStateInitialized_();
        requireVector_(solution, "solution");
        requireVector_(step, "Newton step");

        /*
         * 用户定制 PETSc 在 line-search 前传入正 correction Y；旧 natural
         * limiter 处理的是 delta=-Y。Guard 保证异常路径也会恢复 Y 的符号。
         */
        PetscVectorSignFlipGuard signGuard(
            step);

        const auto &primaryMap =
            backend_.dofMap(
                static_cast<PetscInt>(
                    Indices::numPrimaryVariables));

        const auto &phaseMap =
            backend_.dofMap(
                static_cast<PetscInt>(
                    Indices::numPhaseStateVariables));

        PetscOwnedReadView solutionView(
            solution);

        PetscOwnedWriteView stepView(
            step);

        PetscOwnedReadView phaseView(
            phaseState_);

        for (const Cell cell :
             ownedCells_)
        {
            const PrimaryArray state =
                readOwnedPrimary_(
                    solutionView.data(),
                    solutionView.begin(),
                    solutionView.end(),
                    primaryMap,
                    cell);

            PrimaryArray delta =
                readOwnedPrimary_(
                    stepView.data(),
                    stepView.begin(),
                    stepView.end(),
                    primaryMap,
                    cell);

            // KSP solves in scaled coordinates dU_hat. Recover the physical
            // Natural increment before the existing phase-aware limiter.
            if (options_.scaling.enabled)
            {
                for (int variable = 0; variable < Indices::numPrimaryVariables; ++variable)
                    delta[static_cast<std::size_t>(variable)] *= variableScaleFactor_(variable);
            }

            const auto phaseState =
                readOwnedPhaseState_(
                    phaseView.data(),
                    phaseView.begin(),
                    phaseView.end(),
                    phaseMap,
                    cell);

            if constexpr (Indices::fullyCompositionalThreePhase)
            {
                limitNaturalNewtonIncrement<Indices>(
                    delta, state, phaseState.phasePresence);
                if (options_.maximumAqueousUnsupportedNewtonMoleFraction)
                {
                    limitAqueousCompositionNewtonIncrement<Indices>(
                        delta, state,
                        options_.aqueousNewtonWaterComponent,
                        options_.aqueousNewtonCo2Component,
                        *options_.maximumAqueousUnsupportedNewtonMoleFraction);
                }
            }
            else
            {
                limitNaturalNewtonIncrement<Indices>(
                    delta, state, phaseState.phase);
            }

            writeOwnedPrimary_(
                stepView.data(),
                stepView.begin(),
                stepView.end(),
                primaryMap,
                cell,
                delta);
        }
    }

