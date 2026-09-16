/**
 * @file natural_petsc_runtime_failure_diagnostics.inc
 * @brief Natural/PETSc 非线性失败时的局部质量守恒诊断实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    /**
     * @brief 在失败 Newton 状态上分解 MPI 全局最大组分质量守恒残差。
     *
     * 必须在 rejected attempt 回滚之前调用。函数复用生产 residual 的 cell
     * cache、TPFA face flux 与 well source 公式，对目标单元逐项重算守恒式，
     * 同时记录每个面的 phase potential/upwind/mobility/component flux 以及井
     * 穿孔源项。函数只读求解状态，不修改任何物理/数值变量。
     */
    [[nodiscard]] ResidualFailureDiagnostic evaluateResidualFailureDiagnostic(
        Vec solution,
        Vec residual) const
    {
        requireVector_(solution, "solution");
        requireVector_(residual, "residual");

        ResidualFailureDiagnostic diagnostic;
        diagnostic.timeStep = options_.timeStep;

        const MPI_Comm comm = backend_.communicator();
        PetscMPIInt rank = 0;
        PetscCallMPIAbort(comm, MPI_Comm_rank(comm, &rank));

        const auto &map = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables));
        PetscOwnedReadView residualView(residual);

        // 只在组分质量守恒行中寻找全局最大值；这是本诊断的唯一目标。
        double localMaximum = -1.0;
        Cell localMaximumCell = Cell(-1);
        int localMaximumComponent = -1;
        for (const Cell cell : ownedCells_)
        {
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const int equation = Indices::Equation::massConservation[
                    static_cast<std::size_t>(component)];
                const PetscInt global = map.globalIndex(cell, equation);
                if (global < residualView.begin() || global >= residualView.end())
                    throw std::runtime_error(
                        "Owned Natural mass residual row is outside PETSc ownership range.");
                double value = PetscRealPart(
                    residualView.data()[global - residualView.begin()]);
                if (options_.scaling.enabled)
                    value *= options_.scaling.massResidualScale;
                value = std::abs(value);
                if (value > localMaximum)
                {
                    localMaximum = value;
                    localMaximumCell = cell;
                    localMaximumComponent = component;
                }
            }
        }

        double globalMaximum = -1.0;
        PetscCallMPIAbort(
            comm,
            MPI_Allreduce(&localMaximum, &globalMaximum, 1, MPI_DOUBLE, MPI_MAX, comm));
        if (!(globalMaximum >= 0.0) || !std::isfinite(globalMaximum))
            return diagnostic;

        int candidateRank = localMaximum == globalMaximum
            ? static_cast<int>(rank)
            : std::numeric_limits<int>::max();
        int winningRank = 0;
        PetscCallMPIAbort(
            comm,
            MPI_Allreduce(&candidateRank, &winningRank, 1, MPI_INT, MPI_MIN, comm));

        PetscInt winningCell = -1;
        int winningComponent = -1;
        if (static_cast<int>(rank) == winningRank)
        {
            winningCell = static_cast<PetscInt>(localMaximumCell);
            winningComponent = localMaximumComponent;
        }
        PetscCallMPIAbort(comm, MPI_Bcast(&winningCell, 1, MPIU_INT, winningRank, comm));
        PetscCallMPIAbort(comm, MPI_Bcast(&winningComponent, 1, MPI_INT, winningRank, comm));

        diagnostic.valid = true;
        diagnostic.maximumAbsoluteMassResidual = globalMaximum;
        diagnostic.currentCellId = winningCell;
        diagnostic.componentIndex = winningComponent;

        // 和 production formFunction 完全相同的当前 cache / well source 路径。
        const auto &cache = buildCellCache_(solution);
        const auto &bhp = broadcastWellBottomHolePressures_(cache);
        const auto &wellFunction = evaluateWellFunctionContributions_(cache, bhp);

        if (static_cast<int>(rank) == winningRank)
        {
            const Cell cell = static_cast<Cell>(winningCell);
            const PetscInt block = map.localBlockIndex(cell);
            const auto &entry = cache[static_cast<std::size_t>(block)];
            diagnostic.inputCellId = backend_.inputCellIndex(cell);
            diagnostic.cellVolume = cellVolumeCache_[static_cast<std::size_t>(block)];
            diagnostic.pressure = scalarValue(entry.state.pressure);

            for (int phase = 0; phase < Indices::numPhases; ++phase)
            {
                const std::size_t p = static_cast<std::size_t>(phase);
                diagnostic.saturation[p] = scalarValue(entry.properties.saturation[p]);
                diagnostic.density[p] = scalarValue(entry.properties.density[p]);
                diagnostic.mobility[p] = scalarValue(entry.properties.mobility[p]);
                for (int component = 0; component < Indices::numComponents; ++component)
                {
                    const std::size_t c = static_cast<std::size_t>(component);
                    diagnostic.massFraction[p][c] =
                        scalarValue(entry.properties.massFraction[p][c]);
                }
            }

            const auto currentAccumulation = computeFluidAccumulation<Indices>(
                entry.properties, resolvedCO2Component_());
            const auto &previousProperties =
                previousPropertiesCache_[static_cast<std::size_t>(block)];
            const auto previousAccumulation = computeFluidAccumulation<Indices>(
                previousProperties, resolvedCO2Component_());
            const double storageFactor = diagnostic.cellVolume / options_.timeStep;

            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                diagnostic.currentAccumulationDensity[c] =
                    scalarValue(currentAccumulation.componentMass[c]);
                diagnostic.previousAccumulationDensity[c] =
                    previousAccumulation.componentMass[c];
                diagnostic.accumulationTerm[c] =
                    (diagnostic.currentAccumulationDensity[c] -
                     diagnostic.previousAccumulationDensity[c]) * storageFactor;
            }

            if constexpr (Indices::hasAdsorption)
            {
                const auto currentAdsorption = computeAdsorbedAccumulation<Indices>(
                    entry.properties.porosity,
                    options_.rockDensity,
                    *standardGasDensityPointer_(),
                    entry.properties.adsorbedVolume);
                const auto previousAdsorption = computeAdsorbedAccumulation<Indices>(
                    previousProperties.porosity,
                    options_.rockDensity,
                    *standardGasDensityPointer_(),
                    previousProperties.adsorbedVolume);
                for (int component = 0; component < Indices::numComponents; ++component)
                {
                    const std::size_t c = static_cast<std::size_t>(component);
                    diagnostic.accumulationTerm[c] +=
                        (scalarValue(currentAdsorption[c]) - previousAdsorption[c]) *
                        storageFactor;
                }
            }

            diagnostic.faces.reserve(
                connectionCache_[static_cast<std::size_t>(block)].size());
            for (const auto &connection :
                 connectionCache_[static_cast<std::size_t>(block)])
            {
                NaturalMassBalanceFaceDiagnostic<Indices> face;
                const Cell neighbor = connection.neighborCellId;
                const PetscInt neighborBlock = connection.neighborLocalBlock;
                const auto &exterior = cache[static_cast<std::size_t>(neighborBlock)];

                face.neighborCellId = static_cast<PetscInt>(neighbor);
                face.neighborInputCellId = backend_.inputCellIndex(neighbor);
                face.transmissibility = connection.transmissibility;
                face.gravityTerm = connection.gravityTerm;
                face.neighborPressure = scalarValue(exterior.state.pressure);
                for (int phase = 0; phase < Indices::numPhases; ++phase)
                    face.neighborSaturation[static_cast<std::size_t>(phase)] =
                        scalarValue(exterior.properties.saturation[static_cast<std::size_t>(phase)]);

                for (int phase = 0; phase < Indices::numPhases; ++phase)
                {
                    const std::size_t p = static_cast<std::size_t>(phase);
                    const auto [potential, side] = phasePotentialDifference<Indices, Scalar>(
                        entry.state,
                        entry.properties,
                        exterior.state,
                        exterior.properties,
                        phase,
                        connection.gravityTerm,
                        diagnostic.cellVolume,
                        cellVolumeCache_[static_cast<std::size_t>(neighborBlock)],
                        static_cast<long long>(cell),
                        static_cast<long long>(neighbor));

                    const bool interiorUpwind = side == UpwindSide::Interior;
                    const auto &upwind = interiorUpwind ? entry.properties : exterior.properties;
                    const double mobility = scalarValue(upwind.mobility[p]);
                    const double density = scalarValue(upwind.density[p]);
                    const double potentialValue = scalarValue(potential);
                    const double darcyVolumeRate =
                        -connection.transmissibility * potentialValue * mobility;
                    const double phaseMassRate = darcyVolumeRate * density;

                    face.potentialDifference[p] = potentialValue;
                    face.upwindSide[p] = interiorUpwind ? 0 : 1;
                    face.upwindMobility[p] = mobility;
                    face.upwindDensity[p] = density;
                    face.darcyVolumeRate[p] = darcyVolumeRate;
                    face.phaseMassRate[p] = phaseMassRate;

                    if constexpr (Indices::hasIndependentWaterConservation)
                    {
                        if (phase == Indices::Phase::water)
                        {
                            if constexpr (Indices::hasAqueousCO2Dissolution)
                            {
                                const int co2 = resolvedCO2Component_();
                                const double x = scalarValue(upwind.aqueousCO2MassFraction);
                                face.phaseComponentMassFlux[p][static_cast<std::size_t>(co2)] =
                                    phaseMassRate * x;
                            }
                            continue;
                        }
                    }

                    for (int component = 0; component < Indices::numComponents; ++component)
                    {
                        const std::size_t c = static_cast<std::size_t>(component);
                        face.phaseComponentMassFlux[p][c] =
                            phaseMassRate * scalarValue(upwind.massFraction[p][c]);
                    }
                }

                // 使用生产函数再算一次面总组分通量，避免诊断复制公式造成总量偏差。
                const auto exactFlux = computeFaceMassFlux<Indices, Scalar>(
                    entry.state,
                    entry.properties,
                    exterior.state,
                    exterior.properties,
                    connection.transmissibility,
                    connection.gravityTerm,
                    diagnostic.cellVolume,
                    cellVolumeCache_[static_cast<std::size_t>(neighborBlock)],
                    static_cast<long long>(cell),
                    static_cast<long long>(neighbor),
                    resolvedCO2Component_());
                for (int component = 0; component < Indices::numComponents; ++component)
                {
                    const std::size_t c = static_cast<std::size_t>(component);
                    face.totalComponentMassFlux[c] = scalarValue(exactFlux.component[c]);
                    diagnostic.faceFluxTerm[c] += face.totalComponentMassFlux[c];
                }
                diagnostic.faces.push_back(face);
            }

            const auto &cellWell = wellFunction.cellSource[static_cast<std::size_t>(block)];
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                diagnostic.wellSourceTerm[c] = scalarValue(cellWell.component[c]);
                const int equation = Indices::Equation::massConservation[c];
                const PetscInt global = map.globalIndex(cell, equation);
                diagnostic.petscResidual[c] = PetscRealPart(
                    residualView.data()[global - residualView.begin()]);
                if (options_.scaling.enabled)
                    diagnostic.petscResidual[c] *= options_.scaling.massResidualScale;
                diagnostic.reconstructedResidual[c] =
                    diagnostic.accumulationTerm[c] +
                    diagnostic.faceFluxTerm[c] -
                    diagnostic.wellSourceTerm[c];
                diagnostic.reconstructionError[c] =
                    diagnostic.reconstructedResidual[c] - diagnostic.petscResidual[c];
            }

            for (std::size_t wellIndex = 0; wellIndex < wells_.size(); ++wellIndex)
            {
                const auto &well = wells_[wellIndex];
                if (!well.isActive(currentTime_))
                    continue;
                for (const auto &perforation : well.perforations)
                {
                    if (perforation.currentCellId != cell)
                        continue;

                    NaturalMassBalanceWellDiagnostic<Indices> wellDiagnostic;
                    wellDiagnostic.wellVectorIndex = static_cast<int>(wellIndex);
                    wellDiagnostic.wellId = well.id;
                    wellDiagnostic.wellType = static_cast<int>(well.type);
                    wellDiagnostic.wellControl = static_cast<int>(well.control);
                    wellDiagnostic.active = 1;
                    wellDiagnostic.target = well.target;
                    wellDiagnostic.signedTarget = well.signedTarget();
                    wellDiagnostic.controlledRate = well.isRateControlled()
                        ? wellFunction.controlledRateGlobal[wellIndex]
                        : scalarValue(bhp[wellIndex]);
                    wellDiagnostic.controlResidual = well.isRateControlled()
                        ? wellDiagnostic.controlledRate - wellDiagnostic.signedTarget
                        : scalarValue(bhp[wellIndex]) - well.target;
                    wellDiagnostic.bhp = scalarValue(bhp[wellIndex]);
                    wellDiagnostic.wellIndex = perforation.wellIndex;
                    wellDiagnostic.injectionPhaseFraction = well.injectionPhaseFraction;
                    wellDiagnostic.injectionComponentMassFraction =
                        well.injectionComponentMassFraction;

                    const auto perforationResult = evaluatePerforation_(
                        well, perforation, bhp[wellIndex], entry);
                    for (int phase = 0; phase < Indices::numPhases; ++phase)
                    {
                        const std::size_t p = static_cast<std::size_t>(phase);
                        wellDiagnostic.pressureDrop[p] =
                            scalarValue(entry.state.pressure - bhp[wellIndex]);
                        wellDiagnostic.mobility[p] = scalarValue(entry.properties.mobility[p]);
                        wellDiagnostic.density[p] = scalarValue(entry.properties.density[p]);
                        wellDiagnostic.surfacePhaseRate[p] =
                            scalarValue(perforationResult.surfacePhaseRate[p]);
                        wellDiagnostic.reservoirPhaseRate[p] =
                            scalarValue(perforationResult.reservoirPhaseRate[p]);
                        wellDiagnostic.phaseMassRate[p] =
                            scalarValue(perforationResult.phaseMassRate[p]);
                    }
                    for (int component = 0; component < Indices::numComponents; ++component)
                    {
                        const std::size_t c = static_cast<std::size_t>(component);
                        wellDiagnostic.componentMassSource[c] =
                            scalarValue(perforationResult.componentMassSource[c]);
                    }
                    diagnostic.perforations.push_back(wellDiagnostic);
                }
            }
        }

        // 将 winning rank 的固定数据和可变长 face/perforation 数据广播给 rank 0
        // （也同步到其他 rank，保持 hook 行为一致）。
        PetscCallMPIAbort(comm, MPI_Bcast(&diagnostic.inputCellId, 1, MPIU_INT, winningRank, comm));
        PetscCallMPIAbort(comm, MPI_Bcast(&diagnostic.cellVolume, 1, MPI_DOUBLE, winningRank, comm));
        PetscCallMPIAbort(comm, MPI_Bcast(&diagnostic.pressure, 1, MPI_DOUBLE, winningRank, comm));
        auto bcastDouble = [&](double *data, int count)
        {
            PetscCallMPIAbort(comm, MPI_Bcast(data, count, MPI_DOUBLE, winningRank, comm));
        };
        bcastDouble(diagnostic.saturation.data(), Indices::numPhases);
        bcastDouble(diagnostic.density.data(), Indices::numPhases);
        bcastDouble(diagnostic.mobility.data(), Indices::numPhases);
        for (int phase = 0; phase < Indices::numPhases; ++phase)
            bcastDouble(
                diagnostic.massFraction[static_cast<std::size_t>(phase)].data(),
                Indices::numComponents);
        bcastDouble(diagnostic.currentAccumulationDensity.data(), Indices::numComponents);
        bcastDouble(diagnostic.previousAccumulationDensity.data(), Indices::numComponents);
        bcastDouble(diagnostic.accumulationTerm.data(), Indices::numComponents);
        bcastDouble(diagnostic.faceFluxTerm.data(), Indices::numComponents);
        bcastDouble(diagnostic.wellSourceTerm.data(), Indices::numComponents);
        bcastDouble(diagnostic.reconstructedResidual.data(), Indices::numComponents);
        bcastDouble(diagnostic.petscResidual.data(), Indices::numComponents);
        bcastDouble(diagnostic.reconstructionError.data(), Indices::numComponents);

        const auto broadcastPodVector = [&](auto &values)
        {
            using Value = typename std::decay_t<decltype(values)>::value_type;
            static_assert(std::is_trivially_copyable_v<Value>);
            int count = static_cast<int>(values.size());
            PetscCallMPIAbort(comm, MPI_Bcast(&count, 1, MPI_INT, winningRank, comm));
            if (static_cast<int>(rank) != winningRank)
                values.resize(static_cast<std::size_t>(count));
            const std::size_t byteCount = values.size() * sizeof(Value);
            if (byteCount > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                throw std::runtime_error("Natural failure diagnostic MPI payload is too large.");
            if (byteCount > 0)
            {
                PetscCallMPIAbort(
                    comm,
                    MPI_Bcast(
                        values.data(),
                        static_cast<int>(byteCount),
                        MPI_BYTE,
                        winningRank,
                        comm));
            }
        };
        broadcastPodVector(diagnostic.faces);
        broadcastPodVector(diagnostic.perforations);

        return diagnostic;
    }

