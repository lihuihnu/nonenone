/**
 * @file natural_petsc_runtime_jacobian.inc
 * @brief Natural/PETSc Jacobian 装配及井非局部块实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    PetscErrorCode assembleJacobian_(
        Vec solution,
        Mat matrix)
    {
        PetscFunctionBegin;

        const auto &cache =
            buildCellCache_(solution);

        const auto &bhp =
            broadcastWellBottomHolePressures_(
                cache);

        PetscCall(
            MatZeroEntries(matrix));

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

            /*
             * 所有普通单元都用 wellPressure - value(wellPressure)=0
             * 形成 dummy well-pressure 未知量的单位对角。RATE 井的 BHP
             * 代表单元是唯一例外：该行的 residual 随后被 q-target
             * 完全替换，因此 Jacobian 中原来的 placeholder 单位行也必须
             * 同步替换，而不能与真实 RATE 导数相加。
             */
            if constexpr (Indices::hasWellUnknown)
            {
                if (isRateRepresentative_(localBlock))
                {
                    clearRateWellControlPlaceholder<Indices>(block);
                }
            }

            // local/base block 先写入 self scratch；face transport 的导数随后
            // 直接累加到该 dense block，不再构造/累加 FaceMassFlux<AD> 总和。
            cellJacobianRowBlockScratch_[0] =
                primaryGlobalBlockCache_[static_cast<std::size_t>(localBlock)];
            preparedCellJacobianRowCount_ = 1;
            packCellSelfJacobianBlock_(
                localBlock,
                block);

            prepareFaceJacobianRowsAndAccumulateSelf_(
                cell,
                localBlock,
                cache);

            insertPreparedCellJacobianColumn_(
                matrix,
                localBlock);

        }

        insertWellJacobian_(
            matrix,
            cache,
            bhp);

        PetscCall(
            MatAssemblyBegin(
                matrix,
                MAT_FINAL_ASSEMBLY));

        PetscCall(
            MatAssemblyEnd(
                matrix,
                MAT_FINAL_ASSEMBLY));

        PetscFunctionReturn(PETSC_SUCCESS);
    }

    /**
     * @brief 计算 face AD 导数，并直接写 self/reciprocal Jacobian scratch。
     *
     * Residual 在 G8K 已改为纯标量 face kernel，因此 Jacobian 是每个状态唯一
     * 一遍逐面 AD 传播。这里也不再先构造 `outwardFlux` AD 总和再加到 cell
     * residual；每个 face 的导数直接累加到 self dense block，并同时写 reciprocal
     * row。若同一 CellCache generation 的 residual 已运行，则复用它记录的轻量
     * upwind 决策，保证分支完全一致而不缓存任何逐面 AD 值。
     */
    void prepareFaceJacobianRowsAndAccumulateSelf_(
        Cell cell,
        PetscInt interiorBlock,
        const std::vector<CellCacheEntry> &cache) const
    {
        constexpr PetscInt N =
            static_cast<PetscInt>(Indices::numPrimaryVariables);
        constexpr std::size_t blockValueCount =
            static_cast<std::size_t>(N * N);

        const std::size_t interiorIndex = static_cast<std::size_t>(interiorBlock);
        const auto &connections = connectionCache_[interiorIndex];
        const std::size_t requiredRowCount = connections.size() + std::size_t(1);
        if (cellJacobianRowBlockScratch_.size() < requiredRowCount ||
            cellJacobianBlockValueScratch_.size() < requiredRowCount * blockValueCount)
            throw std::logic_error(
                "Natural cell Jacobian scratch is smaller than the static connection stencil.");

        if (preparedCellJacobianRowCount_ != 1)
            throw std::logic_error(
                "Natural self Jacobian block must be packed before face derivatives.");

        const auto &interior = cache[interiorIndex];
        const double interiorVolume = cellVolumeCache_[interiorIndex];
        const long long interiorCellId = static_cast<long long>(cell);
        const int dissolvedCO2Component = resolvedCO2Component_();

        const bool reuseResidualDecision =
            faceFluxDecisionGeneration_[interiorIndex] == cellCacheGeneration_;
        const auto &decisions = faceFluxDecisionCache_[interiorIndex];
        if (reuseResidualDecision && decisions.size() != connections.size())
            throw std::logic_error(
                "Natural face-flux decision cache does not match static connectivity.");

        PetscScalar *selfValues = cellJacobianBlockValueScratch_.data();

        for (std::size_t connectionIndex = 0;
             connectionIndex < connections.size();
             ++connectionIndex)
        {
            const auto &connection = connections[connectionIndex];
            const Cell neighbor = connection.neighborCellId;
            const PetscInt exteriorBlock = connection.neighborLocalBlock;
            const auto &exterior = cache[static_cast<std::size_t>(exteriorBlock)];

            const auto *decision = reuseResidualDecision
                ? &decisions[connectionIndex]
                : nullptr;

            const auto flux = computeFaceMassFlux<Indices, Scalar>(
                interior.state,
                interior.properties,
                exterior.state,
                exterior.properties,
                connection.transmissibility,
                connection.gravityTerm,
                interiorVolume,
                cellVolumeCache_[static_cast<std::size_t>(exteriorBlock)],
                interiorCellId,
                static_cast<long long>(neighbor),
                dissolvedCO2Component,
                decision);

            const PetscInt scratchRow = preparedCellJacobianRowCount_++;
            cellJacobianRowBlockScratch_[static_cast<std::size_t>(scratchRow)] =
                primaryGlobalBlockCache_[static_cast<std::size_t>(exteriorBlock)];

            PetscScalar *reciprocalValues =
                cellJacobianBlockValueScratch_.data() +
                static_cast<std::size_t>(scratchRow) * blockValueCount;
            std::fill_n(
                reciprocalValues, blockValueCount, PetscScalar(0.0));

            for (int component = 0;
                 component < Indices::numComponents;
                 ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                const PetscInt row = Indices::Equation::massConservation[c];
                const double selfRowScale =
                    equationScaleFactor_(row, interiorBlock);
                const double reciprocalRowScale =
                    equationScaleFactor_(row, exteriorBlock);

                for (PetscInt col = 0; col < N; ++col)
                {
                    const double derivative = flux.component[c].derivative(col);
                    const double variableScale = variableScaleFactor_(col);
                    const std::size_t offset =
                        static_cast<std::size_t>(row * N + col);
                    selfValues[offset] += static_cast<PetscScalar>(
                        derivative * selfRowScale * variableScale);
                    reciprocalValues[offset] = static_cast<PetscScalar>(
                        -derivative * reciprocalRowScale * variableScale);
                }
            }

            if constexpr (Indices::hasIndependentWaterConservation)
            {
                constexpr PetscInt row = Indices::Equation::waterConservation;
                const double selfRowScale =
                    equationScaleFactor_(row, interiorBlock);
                const double reciprocalRowScale =
                    equationScaleFactor_(row, exteriorBlock);
                for (PetscInt col = 0; col < N; ++col)
                {
                    const double derivative = flux.water.derivative(col);
                    const double variableScale = variableScaleFactor_(col);
                    const std::size_t offset =
                        static_cast<std::size_t>(row * N + col);
                    selfValues[offset] += static_cast<PetscScalar>(
                        derivative * selfRowScale * variableScale);
                    reciprocalValues[offset] = static_cast<PetscScalar>(
                        -derivative * reciprocalRowScale * variableScale);
                }
            }
        }
    }

    /** @brief 把 self block 导数写入已准备 scratch 的第 0 个 block。 */
    template <class ResidualBlock>
    void packCellSelfJacobianBlock_(
        PetscInt localBlock,
        const ResidualBlock &block) const
    {
        constexpr PetscInt N =
            static_cast<PetscInt>(Indices::numPrimaryVariables);
        constexpr std::size_t blockValueCount =
            static_cast<std::size_t>(N * N);

        if (preparedCellJacobianRowCount_ < 1 ||
            cellJacobianBlockValueScratch_.size() < blockValueCount)
        {
            throw std::logic_error(
                "Natural cell Jacobian scratch was not prepared before self-block packing.");
        }

        PetscScalar *values = cellJacobianBlockValueScratch_.data();
        for (PetscInt row = 0; row < N; ++row)
        {
            const auto &equation =
                block.value[static_cast<std::size_t>(row)];
            const double rowScale = equationScaleFactor_(row, localBlock);
            for (PetscInt col = 0; col < N; ++col)
            {
                values[static_cast<std::size_t>(row * N + col)] =
                    static_cast<PetscScalar>(
                        equation.derivative(col) *
                        rowScale * variableScaleFactor_(col));
            }
        }
    }

    /**
     * @brief 每个 owned cell 仅调用一次 MatSetValuesBlocked。
     *
     * row list = [self, neighbor_0, neighbor_1, ...]，column list 只有当前
     * interior block，因此不会产生不需要的 row×column 笛卡尔积非零块。
     */
    void insertPreparedCellJacobianColumn_(
        Mat matrix,
        PetscInt interiorBlock) const
    {
        if (preparedCellJacobianRowCount_ <= 0)
            throw std::logic_error("Natural cell Jacobian row scratch is empty.");

        const PetscInt columnBlock =
            primaryGlobalBlockCache_[static_cast<std::size_t>(interiorBlock)];

        PetscCallAbort(
            PetscObjectComm(reinterpret_cast<PetscObject>(matrix)),
            MatSetValuesBlocked(
                matrix,
                preparedCellJacobianRowCount_,
                cellJacobianRowBlockScratch_.data(),
                1,
                &columnBlock,
                cellJacobianBlockValueScratch_.data(),
                ADD_VALUES));
    }

    void insertWellJacobian_(
        [[maybe_unused]] Mat matrix,
        const std::vector<CellCacheEntry> &cache,
        const std::vector<Scalar> &bhp) const
    {
        if constexpr (!Indices::hasWellUnknown)
        {
            return;
        }

        constexpr PetscInt N =
            static_cast<PetscInt>(Indices::numPrimaryVariables);
        constexpr std::size_t rowCapacity =
            static_cast<std::size_t>(Indices::numComponents) +
            (Indices::hasIndependentWaterConservation ? std::size_t(1) : std::size_t(0)) +
            std::size_t(1); // optional RATE control row

        for (std::size_t wellIndex = 0;
             wellIndex < wells_.size();
             ++wellIndex)
        {
            const auto &well = wells_[wellIndex];
            const auto &wellAssembly = wellAssemblyCache_[wellIndex];

            if (wellIndex >= wellBhpGlobalDof_.size() ||
                wellIndex >= wellControlGlobalDof_.size())
            {
                throw std::logic_error(
                    "Natural well global-DOF cache is not initialized.");
            }

            const PetscInt bhpColumn = wellBhpGlobalDof_[wellIndex];
            const PetscInt controlRow = wellControlGlobalDof_[wellIndex];

            if (!well.isActive(currentTime_))
            {
                // Inactive well: no perforation source. Close the additional
                // well-pressure unknown with bhp - p_reservoir = 0.
                const PetscInt localBlock = wellAssembly.representativeLocalBlock;
                if (localBlock >= 0)
                {
                    const PetscInt globalBlock =
                        primaryGlobalBlockCache_[static_cast<std::size_t>(localBlock)];
                    const PetscInt pressureColumn =
                        globalBlock * N + Indices::Primary::pressure;
                    const PetscScalar minusOne = static_cast<PetscScalar>(
                        -wellControlEquationScaleFactor_(wellIndex) *
                        variableScaleFactor_(Indices::Primary::pressure));

                    PetscCallAbort(
                        PetscObjectComm(reinterpret_cast<PetscObject>(matrix)),
                        MatSetValue(
                            matrix,
                            controlRow,
                            pressureColumn,
                            minusOne,
                            ADD_VALUES));
                }
                continue;
            }

            for (const auto &localPerforation : wellAssembly.localPerforations)
            {
                const auto &perforation =
                    well.perforations[localPerforation.perforationIndex];
                const PetscInt localBlock = localPerforation.localBlock;
                const PetscInt globalBlock = localPerforation.globalBlock;

                const auto &result =
                    evaluatePerforation_(
                        well,
                        perforation,
                        bhp[wellIndex],
                        cache[static_cast<std::size_t>(localBlock)]);

                // 所有有效 cols/rows/values 都会在本穿孔调用中完整覆盖。
                // 不做 `{}` 值初始化，避免每个 perforation 先清零整块临时数组。
                std::array<PetscInt, static_cast<std::size_t>(N)> cols;
                const PetscInt firstGlobal = globalBlock * N;
                for (PetscInt col = 0; col < N; ++col)
                {
                    cols[static_cast<std::size_t>(col)] = firstGlobal + col;
                }
                cols[static_cast<std::size_t>(Indices::Primary::wellPressure)] =
                    bhpColumn;

                std::array<PetscInt, rowCapacity> rows;
                std::array<PetscScalar, rowCapacity * static_cast<std::size_t>(N)> values;
                PetscInt rowCount = 0;

                const auto appendScaledRow = [&](
                    PetscInt globalRow,
                    double rowScale,
                    const Scalar &equation)
                {
                    const std::size_t rowOffset =
                        static_cast<std::size_t>(rowCount) * static_cast<std::size_t>(N);
                    rows[static_cast<std::size_t>(rowCount)] = globalRow;
                    for (PetscInt col = 0; col < N; ++col)
                    {
                        values[rowOffset + static_cast<std::size_t>(col)] =
                            static_cast<PetscScalar>(
                                equation.derivative(col) *
                                rowScale * variableScaleFactor_(col));
                    }
                    ++rowCount;
                };

                for (int component = 0;
                     component < Indices::numComponents;
                     ++component)
                {
                    const PetscInt equation =
                        Indices::Equation::massConservation[
                            static_cast<std::size_t>(component)];
                    appendScaledRow(
                        firstGlobal + equation,
                        equationScaleFactor_(equation, localBlock),
                        -result.componentMassSource[static_cast<std::size_t>(component)]);
                }

                if constexpr (Indices::hasIndependentWaterConservation)
                {
                    constexpr PetscInt equation =
                        Indices::Equation::waterConservation;
                    appendScaledRow(
                        firstGlobal + equation,
                        equationScaleFactor_(equation, localBlock),
                        -result.waterMassSource);
                }

                if (well.isRateControlled())
                {
                    const Scalar controlledRate =
                        selectControlledRate<Indices>(
                            well.control,
                            result.surfacePhaseRate,
                            result.reservoirPhaseRate);
                    appendScaledRow(
                        controlRow,
                        wellControlEquationScaleFactor_(wellIndex),
                        controlledRate);
                }

                PetscCallAbort(
                    PetscObjectComm(reinterpret_cast<PetscObject>(matrix)),
                    MatSetValues(
                        matrix,
                        rowCount,
                        rows.data(),
                        N,
                        cols.data(),
                        values.data(),
                        ADD_VALUES));
            }
        }
    }

    void cacheWellGlobalDofs_(
        const std::vector<Well> &wells)
    {
        wellBhpGlobalDof_.assign(wells.size(), PetscInt(-1));
        wellControlGlobalDof_.assign(wells.size(), PetscInt(-1));

        if constexpr (!Indices::hasWellUnknown)
        {
            return;
        }
        else
        {
            if (wells.empty())
            {
                return;
            }

            const auto map = backend_.dofMap(
                static_cast<PetscInt>(Indices::numPrimaryVariables));

            // StructuredGrid builds part of its DOF-layout metadata lazily.
            // That initialization contains PETSc collective operations
            // (for example DMCreateGlobalVector).  Therefore it must be
            // materialized on *all* ranks before the owner-only branch below
            // starts calling map.globalIndex().  Otherwise the owner rank may
            // enter a PETSc collective while the other ranks have already
            // advanced to MPI_Allreduce(), causing an MPI deadlock.
            [[maybe_unused]] const PetscInt globalDofCount =
                map.globalDofCount();

            std::vector<PetscInt> localBhp(wells.size(), PetscInt(-1));
            std::vector<PetscInt> localControl(wells.size(), PetscInt(-1));
            std::vector<PetscInt> ownerCountLocal(wells.size(), PetscInt(0));

            for (std::size_t wellIndex = 0;
                 wellIndex < wells.size();
                 ++wellIndex)
            {
                const Cell bhpCell = wells[wellIndex].bhpCellId;
                if (!backend_.isLocal(bhpCell))
                {
                    continue;
                }

                localBhp[wellIndex] =
                    map.globalIndex(
                        bhpCell,
                        Indices::Primary::wellPressure);

                localControl[wellIndex] =
                    map.globalIndex(
                        bhpCell,
                        Indices::Equation::wellControl);

                ownerCountLocal[wellIndex] = 1;
            }

            std::vector<PetscInt> ownerCountGlobal(
                wells.size(),
                PetscInt(0));

            PetscCallMPIAbort(
                backend_.communicator(),
                MPI_Allreduce(
                    localBhp.data(),
                    wellBhpGlobalDof_.data(),
                    static_cast<int>(wells.size()),
                    MPIU_INT,
                    MPI_MAX,
                    backend_.communicator()));

            PetscCallMPIAbort(
                backend_.communicator(),
                MPI_Allreduce(
                    localControl.data(),
                    wellControlGlobalDof_.data(),
                    static_cast<int>(wells.size()),
                    MPIU_INT,
                    MPI_MAX,
                    backend_.communicator()));

            PetscCallMPIAbort(
                backend_.communicator(),
                MPI_Allreduce(
                    ownerCountLocal.data(),
                    ownerCountGlobal.data(),
                    static_cast<int>(wells.size()),
                    MPIU_INT,
                    MPI_SUM,
                    backend_.communicator()));

            for (std::size_t wellIndex = 0;
                 wellIndex < wells.size();
                 ++wellIndex)
            {
                if (ownerCountGlobal[wellIndex] != 1 ||
                    wellBhpGlobalDof_[wellIndex] < 0 ||
                    wellControlGlobalDof_[wellIndex] < 0)
                {
                    throw std::runtime_error(
                        "Natural well BHP representative cell does not have exactly one owning MPI rank.");
                }
            }
        }
    }

