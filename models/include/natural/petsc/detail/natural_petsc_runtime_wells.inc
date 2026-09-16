/**
 * @file natural_petsc_runtime_wells.inc
 * @brief Natural/PETSc 井 BHP 广播、穿孔源项和井控制残差实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    [[nodiscard]] const std::vector<Scalar> &
    broadcastWellBottomHolePressures_(
        const std::vector<CellCacheEntry> &cache) const
    {
        if (cachedWellBhpCellCacheGeneration_ == cellCacheGeneration_ &&
            wellBhpScalarScratch_.size() == wells_.size())
        {
            return wellBhpScalarScratch_;
        }

        wellBhpLocalScratch_.assign(wells_.size(), 0.0);
        wellBhpGlobalScratch_.assign(wells_.size(), 0.0);

        for (std::size_t wellIndex = 0; wellIndex < wells_.size(); ++wellIndex)
        {
            const PetscInt localBlock =
                wellAssemblyCache_[wellIndex].representativeLocalBlock;
            if (localBlock >= 0)
            {
                wellBhpLocalScratch_[wellIndex] = scalarValue(
                    cache[static_cast<std::size_t>(localBlock)].state.wellPressure);
            }
        }

        if (!wells_.empty())
        {
            PetscCallMPIAbort(
                backend_.communicator(),
                MPI_Allreduce(
                    wellBhpLocalScratch_.data(),
                    wellBhpGlobalScratch_.data(),
                    static_cast<int>(wells_.size()),
                    MPI_DOUBLE, MPI_SUM, backend_.communicator()));
        }

        wellBhpScalarScratch_.assign(wells_.size(), Scalar(0.0));
        for (std::size_t wellIndex = 0; wellIndex < wells_.size(); ++wellIndex)
        {
            wellBhpScalarScratch_[wellIndex] = PrimaryVariables<Indices>::make(
                wellBhpGlobalScratch_[wellIndex],
                Indices::Primary::wellPressure);
        }
        cachedWellBhpCellCacheGeneration_ = cellCacheGeneration_;
        return wellBhpScalarScratch_;
    }

    [[nodiscard]] const PerforationWellResult<
        Indices,
        Scalar> &
    evaluatePerforation_(
        const Well &well,
        const NaturalPerforation<Indices> &perforation,
        const Scalar &bhp,
        const CellCacheEntry &entry) const
    {
        // Runtime assembly is already scratch-based/non-reentrant. Reuse one
        // phase-pressure AD array and one complete perforation result instead of
        // default-constructing both for every perforation in residual/Jacobian.
        for (int phase = 0;
             phase < Indices::numPhases;
             ++phase)
        {
            wellPhasePressureScratch_[
                static_cast<std::size_t>(phase)] =
                entry.state.pressure;
        }

        computePerforationWellSourceInto<Indices, Scalar>(
            wellPerforationResultScratch_,
            wellPerforationWorkspaceScratch_,
            well.type,
            perforation.wellIndex,
            bhp,
            wellPhasePressureScratch_,
            entry.properties.density,
            entry.properties.mobility,
            fluid_.surfaceDensity,
            well.injectionPhaseFraction,
            well.injectionComponentMassFraction,
            entry.properties.massFraction,
            entry.properties.aqueousCO2MassFraction,
            resolvedCO2Component_());
        return wellPerforationResultScratch_;
    }

    [[nodiscard]] const WellFunctionContributions &
    evaluateWellFunctionContributions_(
        const std::vector<CellCacheEntry> &cache,
        const std::vector<Scalar> &bhp) const
    {
        if (wellFunctionScratch_.cellSource.size() != localCacheBlockCount_)
        {
            throw std::logic_error(
                "Natural well-source cache is not initialized for local Vec storage.");
        }

        // 井源项只可能落在穿孔单元。不要在每次 residual 中清零整个 local grid；
        // 只重置静态穿孔 block，使无井/少井大网格的成本与网格规模解耦。
        for (PetscInt localBlock : localWellSourceBlocks_)
        {
            wellFunctionScratch_.cellSource[
                static_cast<std::size_t>(localBlock)] = CellWellSource{};
        }
        wellRateLocalScratch_.assign(wells_.size(), 0.0);
        wellFunctionScratch_.controlledRateGlobal.assign(wells_.size(), 0.0);
        bool hasActiveRateWell = false;

        for (std::size_t wellIndex = 0; wellIndex < wells_.size(); ++wellIndex)
        {
            const auto &well = wells_[wellIndex];
            if (!well.isActive(currentTime_))
                continue;
            hasActiveRateWell = hasActiveRateWell || well.isRateControlled();

            for (const auto &localPerforation :
                 wellAssemblyCache_[wellIndex].localPerforations)
            {
                const auto &perforation =
                    well.perforations[localPerforation.perforationIndex];
                const PetscInt localBlock = localPerforation.localBlock;
                const auto &perforationResult = evaluatePerforation_(
                    well, perforation, bhp[wellIndex],
                    cache[static_cast<std::size_t>(localBlock)]);

                auto &cellSource = wellFunctionScratch_.cellSource[
                    static_cast<std::size_t>(localBlock)];
                for (int component = 0; component < Indices::numComponents; ++component)
                {
                    cellSource.component[static_cast<std::size_t>(component)] +=
                        perforationResult.componentMassSource[
                            static_cast<std::size_t>(component)];
                }
                if constexpr (Indices::hasIndependentWaterConservation)
                    cellSource.water += perforationResult.waterMassSource;

                if (well.isRateControlled())
                {
                    wellRateLocalScratch_[wellIndex] += scalarValue(
                        selectControlledRate<Indices>(
                            well.control, perforationResult.surfacePhaseRate,
                            perforationResult.reservoirPhaseRate));
                }
            }
        }

        if (hasActiveRateWell)
        {
            PetscCallMPIAbort(
                backend_.communicator(),
                MPI_Allreduce(
                    wellRateLocalScratch_.data(),
                    wellFunctionScratch_.controlledRateGlobal.data(),
                    static_cast<int>(wells_.size()),
                    MPI_DOUBLE, MPI_SUM, backend_.communicator()));
        }
        return wellFunctionScratch_;
    }

    void overwriteWellControlResiduals_(
        PetscScalar *residualArray,
        PetscInt residualBegin,
        PetscInt residualEnd,
        const std::vector<Scalar> &bhp,
        const std::vector<double> &controlledRateGlobal,
        const std::vector<CellCacheEntry> &cache) const
    {
        if constexpr (!Indices::hasWellUnknown)
        {
            return;
        }

        for (std::size_t wellIndex = 0;
             wellIndex < wells_.size();
             ++wellIndex)
        {
            const auto &well =
                wells_[wellIndex];

            const PetscInt localBlock =
                wellAssemblyCache_[wellIndex].representativeLocalBlock;
            if (localBlock < 0)
                continue;

            [[maybe_unused]] const PetscInt row =
                wellControlGlobalDof_[wellIndex];

            if (row < residualBegin ||
                row >= residualEnd)
            {
                throw std::runtime_error(
                    "Well-control residual row is not owned by BHP representative rank.");
            }

            const double controlScale = options_.scaling.enabled
                ? equationScaleFactor_(Indices::Equation::wellControl, localBlock)
                : 1.0;

            if (!well.isActive(currentTime_))
            {
                residualArray[
                    row - residualBegin] =
                    static_cast<PetscScalar>(
                        scalarValue(
                            bhp[wellIndex] -
                            cache[
                                static_cast<std::size_t>(
                                    localBlock)]
                                .state.pressure) *
                        controlScale);
            }
            else if (well.isRateControlled())
            {
                residualArray[
                    row - residualBegin] =
                    static_cast<PetscScalar>(
                        controlledRateGlobal[wellIndex] -
                        well.signedTarget()) *
                    static_cast<PetscScalar>(controlScale);
            }
            else
            {
                residualArray[
                    row - residualBegin] =
                    static_cast<PetscScalar>(
                        scalarValue(
                            wellControlResidual(
                                well.control,
                                Scalar(0.0),
                                well.target,
                            bhp[wellIndex])) *
                        controlScale);
            }
        }
    }

    [[nodiscard]] bool isRateRepresentative_(
        PetscInt localBlock) const noexcept
    {
        const std::size_t block = static_cast<std::size_t>(localBlock);
        if (block >= representativeWellByLocalBlock_.size())
            return false;

        const std::size_t wellIndex = representativeWellByLocalBlock_[block];
        if (wellIndex == std::numeric_limits<std::size_t>::max())
            return false;

        const auto &well = wells_[wellIndex];
        return well.isActive(currentTime_) && well.isRateControlled();
    }

    /**
     * @brief 解析井穿孔的 owner/local/global block，供 residual/Jacobian 热路径直接使用。
     *
     * 井几何与穿孔 cell 在 setWells() 后不会随 Newton 改变；控制方式和 schedule
     * 仍从 wells_ 动态读取，因此这里只缓存拓扑/索引，不缓存物理控制状态。
     */
    void cacheWellAssemblyMetadata_(const std::vector<Well> &wells)
    {
        wellAssemblyCache_.assign(wells.size(), WellStaticAssemblyCache{});
        localWellSourceBlocks_.clear();
        representativeWellByLocalBlock_.assign(
            localCacheBlockCount_,
            std::numeric_limits<std::size_t>::max());

        // setup/setWells 路径允许一次性清零整个 local cache；Newton 热循环之后
        // 仅清零 localWellSourceBlocks_ 中实际可能被井修改的条目。
        wellFunctionScratch_.cellSource.assign(
            localCacheBlockCount_, CellWellSource{});

        if (wells.empty())
            return;

        const auto map = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables));
        std::vector<unsigned char> sourceBlockSeen(localCacheBlockCount_, 0);

        for (std::size_t wellIndex = 0; wellIndex < wells.size(); ++wellIndex)
        {
            const auto &well = wells[wellIndex];
            auto &assembly = wellAssemblyCache_[wellIndex];

            if (backend_.isLocal(well.bhpCellId))
            {
                const PetscInt localBlock = map.localBlockIndex(well.bhpCellId);
                assembly.representativeLocalBlock = localBlock;
                representativeWellByLocalBlock_[static_cast<std::size_t>(localBlock)] =
                    wellIndex;
            }

            assembly.localPerforations.reserve(well.perforations.size());
            for (std::size_t perforationIndex = 0;
                 perforationIndex < well.perforations.size();
                 ++perforationIndex)
            {
                const auto &perforation = well.perforations[perforationIndex];
                if (!backend_.isLocal(perforation.currentCellId))
                    continue;

                const PetscInt localBlock =
                    map.localBlockIndex(perforation.currentCellId);
                const PetscInt globalBlock =
                    primaryGlobalBlockCache_[static_cast<std::size_t>(localBlock)];

                assembly.localPerforations.push_back({
                    perforationIndex,
                    localBlock,
                    globalBlock});

                auto &seen = sourceBlockSeen[static_cast<std::size_t>(localBlock)];
                if (seen == 0)
                {
                    seen = 1;
                    localWellSourceBlocks_.push_back(localBlock);
                }
            }
        }
    }

