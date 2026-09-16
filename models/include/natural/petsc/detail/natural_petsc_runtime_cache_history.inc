/**
 * @file natural_petsc_runtime_cache_history.inc
 * @brief Natural/PETSc 运行时静态装配缓存、历史物性和气相历史更新实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    void initializeStaticAssemblyCache_()
    {
        ownedCells_ = backend_.ownedCells();
        localSnapshotCells_ = backend_.localSnapshotCells();

        localCacheBlockCount_ = static_cast<std::size_t>(
            backend_.localCellCount(
                static_cast<PetscInt>(Indices::numPrimaryVariables)));
        if (localCacheBlockCount_ == 0 || localSnapshotCells_.empty())
            throw std::logic_error(
                "Natural local cache does not contain any valid primary blocks.");

        porosityCache_.assign(localCacheBlockCount_, 0.0);
        cellVolumeCache_.assign(localCacheBlockCount_, 0.0);
        maximumGasLocalCache_.assign(localCacheBlockCount_, 0.0);
        connectionCache_.assign(localCacheBlockCount_, {});
        faceFluxDecisionCache_.assign(localCacheBlockCount_, {});
        faceFluxDecisionGeneration_.assign(
            localCacheBlockCount_,
            std::numeric_limits<std::size_t>::max());
        primaryGlobalBlockCache_.assign(localCacheBlockCount_, PetscInt(-1));
        localSnapshotBlockCache_.resize(localSnapshotCells_.size());
        ownedLocalBlockCache_.resize(ownedCells_.size());
        wellFunctionScratch_.cellSource.assign(
            localCacheBlockCount_, CellWellSource{});

        const auto map = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables));
        const auto assembly = backend_.assemblyView();
        std::vector<unsigned char> populated(localCacheBlockCount_, 0);

        for (std::size_t snapshotIndex = 0;
             snapshotIndex < localSnapshotCells_.size();
             ++snapshotIndex)
        {
            const Cell cell = localSnapshotCells_[snapshotIndex];
            const PetscInt mapped = map.localBlockIndex(cell);
            if (mapped < 0 ||
                static_cast<std::size_t>(mapped) >= localCacheBlockCount_)
            {
                throw std::logic_error(
                    "Natural local snapshot contains a block outside local Vec storage.");
            }

            localSnapshotBlockCache_[snapshotIndex] = mapped;

            const std::size_t block = static_cast<std::size_t>(mapped);
            if (populated[block] != 0)
                throw std::logic_error(
                    "Natural local snapshot contains a duplicate local block.");
            populated[block] = 1;

            porosityCache_[block] = assembly.porosity(cell);
            cellVolumeCache_[block] = assembly.cellVolume(cell);

            constexpr PetscInt N =
                static_cast<PetscInt>(Indices::numPrimaryVariables);
            const PetscInt firstGlobal = map.globalIndex(cell, 0);
            if (firstGlobal < 0 || firstGlobal % N != 0)
                throw std::logic_error(
                    "Natural primary DOFs are not aligned to PETSc block rows.");
            for (PetscInt component = 1; component < N; ++component)
            {
                if (map.globalIndex(cell, component) != firstGlobal + component)
                    throw std::logic_error(
                        "Natural primary DOFs are not contiguous within a cell block.");
            }
            primaryGlobalBlockCache_[block] = firstGlobal / N;
        }

        // STAR stencil 的 local Vec 含有未通信的棱/角存储槽；这些槽位可以留空，
        // 但每个 owned 单元以及所有 TPFA 面邻居都必须具有有效缓存条目。
        std::size_t maximumConnectionCount = 0;
        for (std::size_t ownedIndex = 0;
             ownedIndex < ownedCells_.size();
             ++ownedIndex)
        {
            const Cell cell = ownedCells_[ownedIndex];
            const PetscInt block = map.localBlockIndex(cell);
            if (block < 0 ||
                static_cast<std::size_t>(block) >= localCacheBlockCount_ ||
                populated[static_cast<std::size_t>(block)] == 0)
            {
                throw std::logic_error(
                    "Natural owned cell is missing from the communicated local snapshot.");
            }
            ownedLocalBlockCache_[ownedIndex] = block;
        }

        for (std::size_t ownedIndex = 0;
             ownedIndex < ownedCells_.size();
             ++ownedIndex)
        {
            const Cell cell = ownedCells_[ownedIndex];
            const PetscInt block = ownedLocalBlockCache_[ownedIndex];
            const auto connections = assembly.connections(cell);
            auto &cachedConnections =
                connectionCache_[static_cast<std::size_t>(block)];
            cachedConnections.reserve(connections.size());
            maximumConnectionCount =
                std::max(maximumConnectionCount, connections.size());

            for (const auto &connection : connections)
            {
                const PetscInt neighborBlock =
                    map.localBlockIndex(connection.neighborCellId);
                if (neighborBlock < 0 ||
                    static_cast<std::size_t>(neighborBlock) >= localCacheBlockCount_ ||
                    populated[static_cast<std::size_t>(neighborBlock)] == 0)
                {
                    throw std::logic_error(
                        "Natural TPFA face neighbor is missing from the communicated local snapshot.");
                }

                cachedConnections.push_back({
                    connection.neighborCellId,
                    neighborBlock,
                    connection.transmissibility,
                    connection.gravityTerm});
            }

            faceFluxDecisionCache_[static_cast<std::size_t>(block)].resize(
                cachedConnections.size());
        }

        constexpr std::size_t N =
            static_cast<std::size_t>(Indices::numPrimaryVariables);
        const std::size_t maximumJacobianRowCount =
            maximumConnectionCount + std::size_t(1); // self + reciprocal rows
        cellJacobianRowBlockScratch_.resize(maximumJacobianRowCount);
        cellJacobianBlockValueScratch_.resize(
            maximumJacobianRowCount * N * N);
    }

    /**
     * @brief accepted Land history 改变后刷新本地 maximum-gas ghost 快照。
     *
     * maximumGasSaturation_ 只在 commitTimeStep() 更新，因此 Newton trial 内
     * 不需要重复 global-to-local。这个函数是该历史 Vec 的唯一热外 scatter 点。
     */
    void refreshMaximumGasLocalCache_()
    {
        if constexpr (!Indices::hasLandTrapping)
        {
            // Land 关闭时 property evaluator 完全不读取 Sg,max；保持零快照即可。
            return;
        }
        else
        {
            auto maximumGas = backend_.localReadView(1, maximumGasSaturation_);
            ++vecScatterStatistics_.maximumGasSaturation;

            if (maximumGasLocalCache_.size() != localCacheBlockCount_)
                maximumGasLocalCache_.assign(localCacheBlockCount_, 0.0);

            for (PetscInt localBlock : localSnapshotBlockCache_)
            {
                maximumGasLocalCache_[static_cast<std::size_t>(localBlock)] =
                    maximumGas.valueByLocalBlock(localBlock, 0);
            }
        }
    }

    /**
     * @brief 从已接受 history Vec 重算 previous-time properties。
     *
     * maximum-gas 本地快照由 commitTimeStep() 先刷新，因此这里只需同步
     * previous solution 与 previous phase-state 两个真正独立的 ghost Vec。
     */
    void refreshPreviousPropertiesCache_()
    {
        auto previousPrimary = backend_.localReadView(
            static_cast<PetscInt>(Indices::numPrimaryVariables),
            previousSolution_);
        ++vecScatterStatistics_.previousSolution;

        auto previousPhase = backend_.localReadView(
            static_cast<PetscInt>(Indices::numPhaseStateVariables),
            previousPhaseState_);
        ++vecScatterStatistics_.previousPhaseState;

        previousPropertiesCache_.resize(localCacheBlockCount_);
        for (PetscInt localBlock : localSnapshotBlockCache_)
        {
            const std::size_t block = static_cast<std::size_t>(localBlock);
            const auto primaryValues =
                readPrimaryLocal_(previousPrimary, localBlock);
            const auto phaseState =
                readPhaseStateLocal_(previousPhase, localBlock);
            const State state = kernel_.decodeState(primaryValues, phaseState);

            Composition overall{};
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                overall[c] = Scalar(phaseState.overallComposition[c]);
            }

            const auto properties = kernel_.evaluateProperties(
                state,
                Scalar(porosityCache_[block]),
                maximumGasLocalCache_[block],
                options_.landConstant,
                &overall);
            previousPropertiesCache_[block] =
                scalarizeCellProperties<Indices>(properties);
        }

        refreshPreviousAccumulationCache_();
    }

    /**
     * @brief 复用最后一次 accepted residual/Jacobian 的状态快照刷新历史物性。
     *
     * commit 发生时 solution/phase-state 通常正是刚收敛且已构建 CellCache 的
     * Vec。Land maximum 更新只影响物性，不改变 CellState/overall composition，
     * 因此可以避免再次同步 previous solution 与 previous phase-state。
     */
    void refreshPreviousPropertiesFromCurrentCache_()
    {
        if (cellCacheScratch_.size() != localCacheBlockCount_)
            throw std::logic_error(
                "Natural current CellCache cannot seed accepted history properties.");

        previousPropertiesCache_.resize(localCacheBlockCount_);
        for (PetscInt localBlock : localSnapshotBlockCache_)
        {
            const std::size_t block = static_cast<std::size_t>(localBlock);
            const auto &entry = cellCacheScratch_[block];

            if constexpr (!Indices::hasLandTrapping)
            {
                // Land 关闭时 accepted properties 与刚收敛的 current properties
                // 完全相同，直接 scalarize，避免一次额外 EOS/transport 评价。
                previousPropertiesCache_[block] =
                    scalarizeCellProperties<Indices>(entry.properties);
            }
            else
            {
                const auto properties = kernel_.evaluateProperties(
                    entry.state,
                    Scalar(porosityCache_[block]),
                    maximumGasLocalCache_[block],
                    options_.landConstant,
                    &entry.overallComposition);
                previousPropertiesCache_[block] =
                    scalarizeCellProperties<Indices>(properties);
            }
        }

        refreshPreviousAccumulationCache_();
    }

    /**
     * @brief 把 accepted previous properties 压缩成守恒时间项需要的标量缓存。
     *
     * previous state 在整个 Newton solve 内不变；流体/吸附蓄积若在 residual
     * 和 Jacobian 中逐 cell 反复从上一时间层物性重算，是纯重复工作。这里在
     * accepted history 更新时计算一次，热循环只读 component/water 标量。
     */
    void refreshPreviousAccumulationCache_()
    {
        if (previousPropertiesCache_.size() != localCacheBlockCount_)
            throw std::logic_error(
                "Natural previous properties are unavailable for accumulation caching.");

        previousFluidAccumulationCache_.resize(localCacheBlockCount_);
        if constexpr (Indices::hasAdsorption)
            previousAdsorbedAccumulationCache_.resize(localCacheBlockCount_);
        else
            previousAdsorbedAccumulationCache_.clear();

        for (PetscInt localBlock : localSnapshotBlockCache_)
        {
            const std::size_t block = static_cast<std::size_t>(localBlock);
            const auto &properties = previousPropertiesCache_[block];

            previousFluidAccumulationCache_[block] =
                computeFluidAccumulation<Indices>(
                    properties,
                    resolvedCO2Component_());

            if constexpr (Indices::hasAdsorption)
            {
                const auto *standardDensity = standardGasDensityPointer_();
                if (standardDensity == nullptr)
                    throw std::logic_error(
                        "Natural adsorption accumulation cache requires standard gas densities.");

                previousAdsorbedAccumulationCache_[block] =
                    computeAdsorbedAccumulation<Indices>(
                        properties.porosity,
                        options_.rockDensity,
                        *standardDensity,
                        properties.adsorbedVolume);
            }
        }
    }

    [[nodiscard]] bool canReuseCurrentCacheForHistory_(Vec solution) const
    {
        if (!currentCacheValid_ || cachedSolution_ != solution)
            return false;

        PetscObjectState solutionState{};
        PetscObjectState phaseStateState{};
        PetscCallAbort(
            backend_.communicator(),
            VecGetState(solution, &solutionState));
        PetscCallAbort(
            backend_.communicator(),
            VecGetState(phaseState_, &phaseStateState));

        return cachedSolutionState_ == solutionState &&
               cachedPhaseStateState_ == phaseStateState;
    }

    void updateMaximumGasSaturation_(
        Vec solution)
    {
        if constexpr (!Indices::hasLandTrapping)
        {
            (void)solution;
            return;
        }
        else
        {
            const auto &primaryMap =
                backend_.dofMap(
                    static_cast<PetscInt>(
                        Indices::numPrimaryVariables));

            const auto &gasMap =
                backend_.dofMap(1);

            PetscOwnedReadView solutionView(
                solution);

            PetscOwnedWriteView maximumView(
                maximumGasSaturation_);

            for (const Cell cell :
                 ownedCells_)
            {
                const PetscInt vaporGlobal =
                    primaryMap.globalIndex(
                        cell,
                        Indices::Primary::
                            vaporSaturation);

                const PetscInt maxGlobal =
                    gasMap.globalIndex(
                        cell,
                        0);

                if (vaporGlobal < solutionView.begin() ||
                    vaporGlobal >= solutionView.end() ||
                    maxGlobal < maximumView.begin() ||
                    maxGlobal >= maximumView.end())
                {
                    throw std::runtime_error(
                        "Land history Vec ownership is inconsistent with primary Vec.");
                }

                const double vaporSaturation =
                    PetscRealPart(
                        solutionView.data()[
                            vaporGlobal -
                            solutionView.begin()]);

                const PetscInt localIndex =
                    maxGlobal -
                    maximumView.begin();

                const double previousMaximum =
                    PetscRealPart(
                        maximumView.data()[
                            localIndex]);

                maximumView.data()[localIndex] =
                    static_cast<PetscScalar>(
                        std::max(
                            previousMaximum,
                            vaporSaturation));
            }
        }
    }

