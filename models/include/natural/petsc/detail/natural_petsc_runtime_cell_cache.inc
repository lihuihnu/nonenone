/**
 * @file natural_petsc_runtime_cell_cache.inc
 * @brief Natural/PETSc 运行时当前单元物性缓存与面通量累加实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    [[nodiscard]] const std::vector<CellCacheEntry> &
    buildCellCache_(Vec solution) const
    {
        PetscObjectState solutionState{};
        PetscObjectState phaseStateState{};
        PetscObjectState maximumGasState{};
        PetscCallAbort(
            backend_.communicator(),
            VecGetState(
                solution,
                &solutionState));
        PetscCallAbort(
            backend_.communicator(),
            VecGetState(
                phaseState_,
                &phaseStateState));
        if constexpr (Indices::hasLandTrapping)
        {
            PetscCallAbort(
                backend_.communicator(),
                VecGetState(
                    maximumGasSaturation_,
                    &maximumGasState));
        }

        if (currentCacheValid_ &&
            cachedSolution_ == solution &&
            cachedSolutionState_ == solutionState &&
            cachedPhaseStateState_ == phaseStateState &&
            (!Indices::hasLandTrapping ||
             cachedMaximumGasState_ == maximumGasState))
        {
            return cellCacheScratch_;
        }

        auto currentPrimary =
            backend_.localReadView(
                static_cast<PetscInt>(Indices::numPrimaryVariables),
                solution);
        ++vecScatterStatistics_.currentSolution;

        auto currentPhase =
            backend_.localReadView(
                static_cast<PetscInt>(Indices::numPhaseStateVariables),
                phaseState_);
        ++vecScatterStatistics_.currentPhaseState;

        if (previousPropertiesCache_.size() != localCacheBlockCount_)
            throw std::logic_error(
                "Natural previous-property cache is not initialized.");

        cellCacheScratch_.resize(localCacheBlockCount_);
        auto &cache = cellCacheScratch_;

        if (maximumGasLocalCache_.size() != localCacheBlockCount_)
            throw std::logic_error(
                "Natural maximum-gas local cache is not initialized.");

        for (PetscInt localBlock : localSnapshotBlockCache_)
        {
            const std::size_t block = static_cast<std::size_t>(localBlock);

            const auto currentPrimaryValues =
                readPrimaryLocal_(currentPrimary, localBlock);

            const auto currentPhaseState =
                readPhaseStateLocal_(currentPhase, localBlock);

            const State state =
                kernel_.decodeState(
                    currentPrimaryValues,
                    currentPhaseState);

            Composition currentOverall{};

            for (int component = 0;
                 component < Indices::numComponents;
                 ++component)
            {
                const auto c =
                    static_cast<std::size_t>(component);

                currentOverall[c] =
                    Scalar(
                        currentPhaseState
                            .overallComposition[c]);

            }

            const double maxGas =
                maximumGasLocalCache_[block];

            const double porosity =
                porosityCache_[block];

            const Properties currentProperties =
                kernel_.evaluateProperties(
                    state,
                    Scalar(porosity),
                    maxGas,
                    options_.landConstant,
                    &currentOverall);

            auto &entry = cache[block];

            entry.state = state;
            entry.overallComposition = currentOverall;
            entry.properties = currentProperties;
        }

        cachedSolution_ = solution;
        cachedSolutionState_ = solutionState;
        cachedPhaseStateState_ = phaseStateState;
        if constexpr (Indices::hasLandTrapping)
            cachedMaximumGasState_ = maximumGasState;
        ++cellCacheGeneration_;
        currentCacheValid_ = true;
        return cache;
    }

    void accumulateOutwardFaceFluxValue_(
        Cell cell,
        PetscInt interiorBlock,
        const std::vector<CellCacheEntry> &cache,
        FaceMassFlux<Indices, double> &sum) const
    {
        const std::size_t blockIndex =
            static_cast<std::size_t>(interiorBlock);
        const auto &interior = cache[blockIndex];
        const double interiorVolume = cellVolumeCache_[blockIndex];
        const long long interiorCellId = static_cast<long long>(cell);
        const int dissolvedCO2Component = resolvedCO2Component_();

        const auto &connections = connectionCache_[blockIndex];
        auto &decisions = faceFluxDecisionCache_[blockIndex];
        if (decisions.size() != connections.size())
            throw std::logic_error(
                "Natural face-flux decision cache does not match static connectivity.");

        for (std::size_t connectionIndex = 0;
             connectionIndex < connections.size();
             ++connectionIndex)
        {
            const auto &connection = connections[connectionIndex];
            const Cell neighbor = connection.neighborCellId;
            const PetscInt exteriorBlock = connection.neighborLocalBlock;
            const auto &exterior = cache[static_cast<std::size_t>(exteriorBlock)];

            const auto flux = computeFaceMassFluxValue<Indices, Scalar>(
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
                &decisions[connectionIndex]);

            for (int component = 0;
                 component < Indices::numComponents;
                 ++component)
            {
                sum.component[static_cast<std::size_t>(component)] +=
                    flux.component[static_cast<std::size_t>(component)];
            }

            if constexpr (Indices::hasIndependentWaterConservation)
                sum.water += flux.water;
        }

        // 只有完整扫完该 owned cell 的所有连接后才发布 generation，避免异常路径
        // 留下半填充的离散分支供 Jacobian 误用。
        faceFluxDecisionGeneration_[blockIndex] = cellCacheGeneration_;
    }

