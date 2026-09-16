/**
 * @file natural_petsc_runtime_residual.inc
 * @brief Natural/PETSc SNES 残差与 Jacobian 公共入口实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    /** @brief 兼容直接调用的 residual 入口。 */
    PetscErrorCode formFunction(
        Vec solution,
        Vec residual)
    {
        return formFunction(nullptr, solution, residual);
    }

    /**
     * @brief 标准 SNES residual 实现；相态更新失败时向 SNES 报 function-domain error。
     */
    PetscErrorCode formFunction(
        SNES snes,
        Vec solution,
        Vec residual)
    {
        PetscFunctionBegin;

        requireReadyForResidual_();
        requireVector_(solution, "solution");
        requireVector_(residual, "residual");

        if (pendingThermoStateFailure_)
        {
            PetscCall(VecSet(residual, 0.0));
            if (snes != nullptr)
            {
                PetscCall(SNESSetFunctionDomainError(snes));
                PetscFunctionReturn(PETSC_SUCCESS);
            }
            SETERRQ(
                backend_.communicator(),
                PETSC_ERR_ARG_OUTOFRANGE,
                "Natural phase-state update reported a recoverable thermodynamic failure before residual evaluation.");
        }

        const auto &cache =
            buildCellCache_(solution);

        const auto &bhp =
            broadcastWellBottomHolePressures_(
                cache);

        const auto &wellFunction =
            evaluateWellFunctionContributions_(
                cache,
                bhp);

        PetscOwnedWriteView residualView(
            residual);

        for (std::size_t ownedIndex = 0;
             ownedIndex < ownedCells_.size();
             ++ownedIndex)
        {
            const Cell cell = ownedCells_[ownedIndex];
            const PetscInt localBlock =
                ownedLocalBlockCache_[ownedIndex];

            const auto &entry =
                cache[
                    static_cast<std::size_t>(
                        localBlock)];

            FaceMassFlux<Indices, double>
                outwardFlux{};

            accumulateOutwardFaceFluxValue_(
                cell,
                localBlock,
                cache,
                outwardFlux);

            const auto &cellWell =
                wellFunction.cellSource[
                    static_cast<std::size_t>(
                        localBlock)];

            const std::size_t blockIndex =
                static_cast<std::size_t>(localBlock);
            const auto *previousAdsorption =
                previousAdsorbedAccumulationCache_.empty()
                    ? nullptr
                    : &previousAdsorbedAccumulationCache_[blockIndex];

            auto block =
                kernel_.assembleLocalResidual(
                    entry.state,
                    entry.properties,
                    previousPropertiesCache_[blockIndex],
                    cellVolumeCache_[blockIndex],
                    options_.timeStep,
                    options_.fugacityScalingFactor,
                    resolvedCO2Component_(),
                    options_.rockDensity,
                    standardGasDensityPointer_(),
                    &previousFluidAccumulationCache_[blockIndex],
                    previousAdsorption);

            addCellConservationTransportAndWell<Indices>(
                block,
                outwardFlux,
                cellWell.component,
                cellWell.water);

            writeOwnedResidualBlock_(
                residualView.data(),
                residualView.begin(),
                residualView.end(),
                primaryGlobalBlockCache_[static_cast<std::size_t>(localBlock)],
                localBlock,
                block);
        }

        overwriteWellControlResiduals_(
            residualView.data(),
            residualView.begin(),
            residualView.end(),
            bhp,
            wellFunction.controlledRateGlobal,
            cache);

        residualView.restore();

        PetscFunctionReturn(PETSC_SUCCESS);
    }

    /** 标准 SNES Jacobian 实现。 */
    PetscErrorCode formJacobian(
        Vec solution,
        Mat jacobian,
        Mat preconditioner)
    {
        PetscFunctionBegin;

        static_assert(
            Indices::useAutomaticDifferentiation,
            "Natural PETSc analytic Jacobian requires ADIndices.");

        requireReadyForResidual_();
        requireVector_(solution, "solution");

        PetscCheck(
            jacobian != nullptr,
            backend_.communicator(),
            PETSC_ERR_ARG_NULL,
            "Jacobian matrix must not be null.");

        PetscCall(
            assembleJacobian_(
                solution,
                jacobian));

        // G8I 在 block/scalar 打包阶段直接应用 row/column scaling。
        // 因此这里不再刷新两个 scaling Vec，也不再整矩阵 MatDiagonalScale。
        if (preconditioner != nullptr &&
            preconditioner != jacobian)
        {
            PetscCall(
                assembleJacobian_(
                    solution,
                    preconditioner));
        }

        PetscFunctionReturn(PETSC_SUCCESS);
    }

