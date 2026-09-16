/**
 * @file natural_petsc_runtime_output_diagnostics.inc
 * @brief Natural/PETSc 全局井状态、库存和储层输出诊断实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    /**
     * @brief 计算 MPI 全局一致的收敛井状态供输出与控制使用。
     *
     * The returned vector is aligned with wells(). Phase rates preserve the
     * simulator sign convention: injection is positive and production is
     * negative. bottomHolePressure and all phase rates are MPI-global values.
     */
    [[nodiscard]] std::vector<WellStateType> evaluateWellStates(
        Vec solution) const
    {
        requireReadyForResidual_();
        requireVector_(solution, "solution");

        std::vector<WellStateType> states(wells_.size());
        if (wells_.empty())
            return states;

        const auto &cache = buildCellCache_(solution);
        const auto &bhp = broadcastWellBottomHolePressures_(cache);

        constexpr std::size_t phaseCount =
            static_cast<std::size_t>(Indices::numPhases);
        constexpr std::size_t componentCount =
            static_cast<std::size_t>(Indices::numComponents);
        constexpr std::size_t componentOffset = 3 * phaseCount;
        constexpr std::size_t valuesPerWell = componentOffset + componentCount;

        std::vector<double> local(wells_.size() * valuesPerWell, 0.0);
        std::vector<double> global(wells_.size() * valuesPerWell, 0.0);

        const auto &map = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables));

        for (std::size_t wellIndex = 0; wellIndex < wells_.size(); ++wellIndex)
        {
            const auto &well = wells_[wellIndex];
            const std::size_t base = wellIndex * valuesPerWell;

            if (!well.isActive(currentTime_))
                continue;

            for (const auto &perforation : well.perforations)
            {
                const Cell cell = perforation.currentCellId;
                if (!backend_.isLocal(cell))
                    continue;

                const PetscInt localBlock = map.localBlockIndex(cell);
                const auto result = evaluatePerforation_(
                    well,
                    perforation,
                    bhp[wellIndex],
                    cache[static_cast<std::size_t>(localBlock)]);

                for (std::size_t phase = 0; phase < phaseCount; ++phase)
                {
                    local[base + phase] +=
                        scalarValue(result.surfacePhaseRate[phase]);
                    local[base + phaseCount + phase] +=
                        scalarValue(result.reservoirPhaseRate[phase]);
                    local[base + 2 * phaseCount + phase] +=
                        scalarValue(result.phaseMassRate[phase]);
                }
                for (std::size_t component = 0;
                     component < componentCount; ++component)
                {
                    local[base + componentOffset + component] +=
                        scalarValue(result.componentMassSource[component]);
                }
            }
        }

        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(
                local.data(),
                global.data(),
                static_cast<int>(global.size()),
                MPI_DOUBLE,
                MPI_SUM,
                backend_.communicator()));

        for (std::size_t wellIndex = 0; wellIndex < wells_.size(); ++wellIndex)
        {
            auto &state = states[wellIndex];
            state.bottomHolePressure = scalarValue(bhp[wellIndex]);
            const std::size_t base = wellIndex * valuesPerWell;

            for (std::size_t phase = 0; phase < phaseCount; ++phase)
            {
                state.surfacePhaseRate[phase] = global[base + phase];
                state.reservoirPhaseRate[phase] =
                    global[base + phaseCount + phase];
                state.phaseMassRate[phase] =
                    global[base + 2 * phaseCount + phase];

                const double reservoirRate = state.reservoirPhaseRate[phase];
                state.flowWeightedDensity[phase] =
                    std::abs(reservoirRate) > 1.0e-30
                        ? state.phaseMassRate[phase] / reservoirRate
                        : 0.0;
            }
            for (std::size_t component = 0;
                 component < componentCount; ++component)
            {
                state.componentMassRate[component] =
                    global[base + componentOffset + component];
            }
        }

        return states;
    }

    /**
     * @brief 计算 MPI 全局守恒组分井源汇率 [kg/s]。
     *
     * The calculation reuses the exact perforation source kernel used by the
     * residual.  Positive contributions are accumulated as injection and
     * negative contributions as positive production magnitudes.  Consequently
     * `netSource = injected - produced` component by component.  This is the
     * quantity required for discrete global mass-balance accounting.
     */
    [[nodiscard]] NaturalGlobalWellComponentRates<Indices>
    evaluateGlobalWellComponentRates(Vec solution) const
    {
        requireReadyForResidual_();
        requireVector_(solution, "solution");

        constexpr std::size_t n =
            static_cast<std::size_t>(Indices::numComponents);
        NaturalGlobalWellComponentRates<Indices> result;
        if (wells_.empty())
            return result;

        const auto &cache = buildCellCache_(solution);
        const auto &bhp = broadcastWellBottomHolePressures_(cache);
        const auto &map = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables));

        constexpr std::size_t offsetNet = 0;
        constexpr std::size_t offsetInjected = offsetNet + n;
        constexpr std::size_t offsetProduced = offsetInjected + n;
        constexpr std::size_t offsetWaterNet = offsetProduced + n;
        constexpr std::size_t offsetWaterInjected = offsetWaterNet + 1;
        constexpr std::size_t offsetWaterProduced = offsetWaterInjected + 1;
        constexpr std::size_t count = offsetWaterProduced + 1;
        std::array<double, count> local{};
        std::array<double, count> global{};

        for (std::size_t wellIndex = 0; wellIndex < wells_.size(); ++wellIndex)
        {
            const auto &well = wells_[wellIndex];
            if (!well.isActive(currentTime_))
                continue;

            for (const auto &perforation : well.perforations)
            {
                const Cell cell = perforation.currentCellId;
                if (!backend_.isLocal(cell))
                    continue;

                const PetscInt localBlock = map.localBlockIndex(cell);
                const auto perforationResult = evaluatePerforation_(
                    well,
                    perforation,
                    bhp[wellIndex],
                    cache[static_cast<std::size_t>(localBlock)]);

                for (std::size_t c = 0; c < n; ++c)
                {
                    const double source =
                        scalarValue(perforationResult.componentMassSource[c]);
                    local[offsetNet + c] += source;
                    if (source >= 0.0)
                        local[offsetInjected + c] += source;
                    else
                        local[offsetProduced + c] += -source;
                }

                if constexpr (Indices::hasIndependentWaterConservation)
                {
                    const double source =
                        scalarValue(perforationResult.waterMassSource);
                    local[offsetWaterNet] += source;
                    if (source >= 0.0)
                        local[offsetWaterInjected] += source;
                    else
                        local[offsetWaterProduced] += -source;
                }
            }
        }

        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(
                local.data(),
                global.data(),
                static_cast<int>(count),
                MPI_DOUBLE,
                MPI_SUM,
                backend_.communicator()));

        for (std::size_t c = 0; c < n; ++c)
        {
            result.netSource[c] = global[offsetNet + c];
            result.injected[c] = global[offsetInjected + c];
            result.produced[c] = global[offsetProduced + c];
        }
        if constexpr (Indices::hasIndependentWaterConservation)
        {
            result.waterNetSource = global[offsetWaterNet];
            result.waterInjected = global[offsetWaterInjected];
            result.waterProduced = global[offsetWaterProduced];
        }
        return result;
    }

    /**
     * @brief 计算当前已收敛状态的 MPI 全局质量库存。
     *
     * 这项统计与网格类型无关，可用于 StructuredGrid/CpGrid：
     * - fluidComponentMass：油/气/水中非水组分总质量（溶解 CO2 已包含）；
     * - waterMass：水相中的 H2O 质量；
     * - dissolvedCO2Mass：水相 CO2 的独立统计；
     * - trappedGasMass / trappedComponentMass：Land 残余气库存；
     * - adsorbedComponentMass：竞争吸附库存。
     */
    [[nodiscard]] NaturalGlobalInventory<Indices>
    evaluateGlobalInventory(Vec solution) const
    {
        requireReadyForResidual_();
        requireVector_(solution, "solution");

        const auto &cache = buildCellCache_(solution);
        const auto &map = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables));

        constexpr std::size_t n =
            static_cast<std::size_t>(Indices::numComponents);
        constexpr std::size_t offsetFluid = 0;
        constexpr std::size_t offsetOilPhase = offsetFluid + n;
        constexpr std::size_t offsetWater = offsetOilPhase + n;
        constexpr std::size_t offsetDissolved = offsetWater + 1;
        constexpr std::size_t offsetTrappedTotal = offsetDissolved + 1;
        constexpr std::size_t offsetTrappedComponents = offsetTrappedTotal + 1;
        constexpr std::size_t offsetAdsorbed = offsetTrappedComponents + n;
        constexpr std::size_t count = offsetAdsorbed + n;

        std::array<double, count> local{};
        std::array<double, count> global{};

        for (const Cell cell : ownedCells_)
        {
            const PetscInt localBlock = map.localBlockIndex(cell);
            const auto &properties =
                cache[static_cast<std::size_t>(localBlock)].properties;
            const double volume =
                cellVolumeCache_[static_cast<std::size_t>(localBlock)];

            const auto accumulation =
                computeFluidAccumulation<Indices>(
                    properties,
                    resolvedCO2Component_());

            for (int component = 0;
                 component < Indices::numComponents;
                 ++component)
            {
                const std::size_t c =
                    static_cast<std::size_t>(component);
                local[offsetFluid + c] +=
                    scalarValue(accumulation.componentMass[c]) * volume;
            }

            {
                const std::size_t liquid =
                    static_cast<std::size_t>(Indices::Phase::liquid);

                const double oilMass =
                    scalarValue(
                        properties.porosity *
                        properties.density[liquid] *
                        properties.saturation[liquid]) *
                    volume;

                for (int component = 0;
                     component < Indices::numComponents;
                     ++component)
                {
                    const std::size_t c =
                        static_cast<std::size_t>(component);

                    local[offsetOilPhase + c] +=
                        oilMass *
                        scalarValue(
                            properties.massFraction[liquid][c]);
                }
            }

            if constexpr (Indices::hasIndependentWaterConservation)
            {
                local[offsetWater] +=
                    scalarValue(accumulation.waterMass) * volume;
            }

            if constexpr (Indices::hasAqueousCO2Dissolution)
            {
                const std::size_t water =
                    static_cast<std::size_t>(Indices::Phase::water);
                const double aqueousCO2Mass =
                    scalarValue(
                        properties.porosity *
                        properties.density[water] *
                        properties.saturation[water] *
                        properties.aqueousCO2MassFraction) *
                    volume;
                local[offsetDissolved] += aqueousCO2Mass;
            }

            if constexpr (Indices::hasLandTrapping)
            {
                const std::size_t vapor =
                    static_cast<std::size_t>(Indices::Phase::vapor);

                const double trappedMass =
                    scalarValue(
                        properties.porosity *
                        properties.density[vapor] *
                        properties.trappedGasSaturation) *
                    volume;

                local[offsetTrappedTotal] += trappedMass;

                for (int component = 0;
                     component < Indices::numComponents;
                     ++component)
                {
                    const std::size_t c =
                        static_cast<std::size_t>(component);
                    local[offsetTrappedComponents + c] +=
                        trappedMass *
                        scalarValue(
                            properties.massFraction[vapor][c]);
                }
            }

            if constexpr (Indices::hasAdsorption)
            {
                const auto *rhoStd = standardGasDensityPointer_();
                if (rhoStd == nullptr)
                    throw std::logic_error(
                        "Adsorption inventory requires standard gas density.");

                const auto adsorbed =
                    computeAdsorbedAccumulation<Indices>(
                        properties.porosity,
                        options_.rockDensity,
                        *rhoStd,
                        properties.adsorbedVolume);

                for (int component = 0;
                     component < Indices::numComponents;
                     ++component)
                {
                    const std::size_t c =
                        static_cast<std::size_t>(component);
                    local[offsetAdsorbed + c] +=
                        scalarValue(adsorbed[c]) * volume;
                }
            }
        }

        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(
                local.data(),
                global.data(),
                static_cast<int>(count),
                MPI_DOUBLE,
                MPI_SUM,
                backend_.communicator()));

        NaturalGlobalInventory<Indices> result;

        for (int component = 0;
             component < Indices::numComponents;
             ++component)
        {
            const std::size_t c =
                static_cast<std::size_t>(component);
            result.fluidComponentMass[c] =
                global[offsetFluid + c];
            result.oilPhaseComponentMass[c] =
                global[offsetOilPhase + c];
            result.trappedComponentMass[c] =
                global[offsetTrappedComponents + c];
            result.adsorbedComponentMass[c] =
                global[offsetAdsorbed + c];
        }

        result.waterMass = global[offsetWater];
        result.dissolvedCO2Mass = global[offsetDissolved];
        result.trappedGasMass = global[offsetTrappedTotal];
        return result;
    }


    /**
     * @brief 计算用于输出/调试的紧凑 MPI 全局储层诊断量。
     */
    [[nodiscard]] NaturalGlobalDiagnostics<Indices>
    evaluateGlobalDiagnostics(Vec solution) const
    {
        requireReadyForResidual_();
        requireVector_(solution, "solution");

        const auto &cache = buildCellCache_(solution);
        const auto &map = backend_.dofMap(
            static_cast<PetscInt>(Indices::numPrimaryVariables));

        constexpr std::size_t phaseCount =
            static_cast<std::size_t>(Indices::numPhases);
        constexpr double activeSaturation = 1.0e-12;
        const double inf = std::numeric_limits<double>::infinity();

        long long localCellCount = 0;
        double localBulkVolume = 0.0;
        double localPoreVolume = 0.0;
        double localPressureSum = 0.0;
        double localPorositySum = 0.0;
        double localPressureMin = inf;
        double localPressureMax = -inf;
        double localPorosityMin = inf;
        double localPorosityMax = -inf;

        std::array<double, phaseCount> localPhaseVolume{};
        std::array<double, phaseCount> localSaturationSum{};
        std::array<double, phaseCount> localSaturationMin{};
        std::array<double, phaseCount> localSaturationMax{};
        std::array<double, phaseCount> localDensitySum{};
        std::array<double, phaseCount> localDensityMin{};
        std::array<double, phaseCount> localDensityMax{};
        std::array<double, phaseCount> localViscositySum{};
        std::array<double, phaseCount> localViscosityMin{};
        std::array<double, phaseCount> localViscosityMax{};

        localSaturationMin.fill(inf);
        localSaturationMax.fill(-inf);
        localDensityMin.fill(inf);
        localDensityMax.fill(-inf);
        localViscosityMin.fill(inf);
        localViscosityMax.fill(-inf);

        long long localLiquidOnly = 0;
        long long localVaporOnly = 0;
        long long localTwoPhase = 0;
        std::array<long long, 8> localPhasePresence{};
        std::array<long long, 3> localPhaseSuppression{};
        double localAqueousWeight = 0.0;
        double localAqueousCO2Sum = 0.0;
        double localAqueousCO2Max = 0.0;
        double localTrappedSum = 0.0;
        double localTrappedMax = 0.0;

        for (const Cell cell : ownedCells_)
        {
            ++localCellCount;
            const PetscInt localBlock = map.localBlockIndex(cell);
            const auto &entry = cache[static_cast<std::size_t>(localBlock)];
            const auto &state = entry.state;
            const auto &properties = entry.properties;

            const double volume =
                cellVolumeCache_[static_cast<std::size_t>(localBlock)];
            const double porosity = scalarValue(properties.porosity);
            const double poreVolume = porosity * volume;
            const double pressure = scalarValue(state.pressure);

            localBulkVolume += volume;
            localPoreVolume += poreVolume;
            localPressureSum += pressure * volume;
            localPorositySum += porosity * volume;
            localPressureMin = std::min(localPressureMin, pressure);
            localPressureMax = std::max(localPressureMax, pressure);
            localPorosityMin = std::min(localPorosityMin, porosity);
            localPorosityMax = std::max(localPorosityMax, porosity);

            for (std::size_t phase = 0; phase < phaseCount; ++phase)
            {
                const double saturation = scalarValue(properties.saturation[phase]);
                const double phaseVolume = std::max(0.0, saturation) * poreVolume;
                localSaturationSum[phase] += saturation * poreVolume;
                localSaturationMin[phase] = std::min(localSaturationMin[phase], saturation);
                localSaturationMax[phase] = std::max(localSaturationMax[phase], saturation);
                localPhaseVolume[phase] += phaseVolume;

                if (saturation > activeSaturation)
                {
                    const double density = scalarValue(properties.density[phase]);
                    const double viscosity = scalarValue(properties.viscosity[phase]);
                    localDensitySum[phase] += density * phaseVolume;
                    localViscositySum[phase] += viscosity * phaseVolume;
                    localDensityMin[phase] = std::min(localDensityMin[phase], density);
                    localDensityMax[phase] = std::max(localDensityMax[phase], density);
                    localViscosityMin[phase] = std::min(localViscosityMin[phase], viscosity);
                    localViscosityMax[phase] = std::max(localViscosityMax[phase], viscosity);
                }
            }

            if constexpr (Indices::fullyCompositionalThreePhase)
            {
                const auto mask = static_cast<std::size_t>(state.phasePresence.bits());
                if (mask < localPhasePresence.size())
                    ++localPhasePresence[mask];

                const auto suppression = state.phaseSuppression;
                if (suppression.contains(CompositionalPhase::Oil))
                    ++localPhaseSuppression[0];
                if (suppression.contains(CompositionalPhase::Gas))
                    ++localPhaseSuppression[1];
                if (suppression.contains(CompositionalPhase::Water))
                    ++localPhaseSuppression[2];
            }
            else
            {
                switch (state.hydrocarbonPhaseState)
                {
                case HydrocarbonPhaseState::LiquidOnly: ++localLiquidOnly; break;
                case HydrocarbonPhaseState::VaporOnly:  ++localVaporOnly; break;
                case HydrocarbonPhaseState::TwoPhase:   ++localTwoPhase; break;
                }
            }

            if constexpr (Indices::hasAqueousCO2Dissolution)
            {
                const std::size_t water = static_cast<std::size_t>(Indices::Phase::water);
                const double waterVolume =
                    std::max(0.0, scalarValue(properties.saturation[water])) * poreVolume;
                const double massFraction = scalarValue(properties.aqueousCO2MassFraction);
                localAqueousWeight += waterVolume;
                localAqueousCO2Sum += massFraction * waterVolume;
                localAqueousCO2Max = std::max(localAqueousCO2Max, massFraction);
            }

            if constexpr (Indices::hasLandTrapping)
            {
                const double trapped = scalarValue(properties.trappedGasSaturation);
                localTrappedSum += trapped * poreVolume;
                localTrappedMax = std::max(localTrappedMax, trapped);
            }
        }

        NaturalGlobalDiagnostics<Indices> result;
        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(&localCellCount, &result.cellCount, 1, MPI_LONG_LONG, MPI_SUM,
                          backend_.communicator()));

        std::array<double, 7> localScalarSum{
            localBulkVolume, localPoreVolume, localPressureSum, localPorositySum,
            localAqueousWeight, localAqueousCO2Sum, localTrappedSum};
        std::array<double, 7> globalScalarSum{};
        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(localScalarSum.data(), globalScalarSum.data(), 7, MPI_DOUBLE, MPI_SUM,
                          backend_.communicator()));

        std::array<double, 2> localMin{localPressureMin, localPorosityMin};
        std::array<double, 2> globalMin{};
        std::array<double, 4> localMax{
            localPressureMax, localPorosityMax, localAqueousCO2Max, localTrappedMax};
        std::array<double, 4> globalMax{};
        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(localMin.data(), globalMin.data(), 2, MPI_DOUBLE, MPI_MIN,
                          backend_.communicator()));
        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(localMax.data(), globalMax.data(), 4, MPI_DOUBLE, MPI_MAX,
                          backend_.communicator()));

        auto reducePhase = [&](const auto &local, auto &global, MPI_Op operation) {
            PetscCallMPIAbort(
                backend_.communicator(),
                MPI_Allreduce(local.data(), global.data(), static_cast<int>(phaseCount),
                              MPI_DOUBLE, operation, backend_.communicator()));
        };

        std::array<double, phaseCount> globalPhaseVolume{};
        std::array<double, phaseCount> globalSaturationSum{};
        std::array<double, phaseCount> globalSaturationMin{};
        std::array<double, phaseCount> globalSaturationMax{};
        std::array<double, phaseCount> globalDensitySum{};
        std::array<double, phaseCount> globalDensityMin{};
        std::array<double, phaseCount> globalDensityMax{};
        std::array<double, phaseCount> globalViscositySum{};
        std::array<double, phaseCount> globalViscosityMin{};
        std::array<double, phaseCount> globalViscosityMax{};

        reducePhase(localPhaseVolume, globalPhaseVolume, MPI_SUM);
        reducePhase(localSaturationSum, globalSaturationSum, MPI_SUM);
        reducePhase(localSaturationMin, globalSaturationMin, MPI_MIN);
        reducePhase(localSaturationMax, globalSaturationMax, MPI_MAX);
        reducePhase(localDensitySum, globalDensitySum, MPI_SUM);
        reducePhase(localDensityMin, globalDensityMin, MPI_MIN);
        reducePhase(localDensityMax, globalDensityMax, MPI_MAX);
        reducePhase(localViscositySum, globalViscositySum, MPI_SUM);
        reducePhase(localViscosityMin, globalViscosityMin, MPI_MIN);
        reducePhase(localViscosityMax, globalViscosityMax, MPI_MAX);

        std::array<long long, 3> localStates{localLiquidOnly, localVaporOnly, localTwoPhase};
        std::array<long long, 3> globalStates{};
        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(localStates.data(), globalStates.data(), 3, MPI_LONG_LONG, MPI_SUM,
                          backend_.communicator()));

        std::array<long long, 8> globalPhasePresence{};
        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(localPhasePresence.data(), globalPhasePresence.data(),
                          static_cast<int>(localPhasePresence.size()), MPI_LONG_LONG, MPI_SUM,
                          backend_.communicator()));

        std::array<long long, 3> globalPhaseSuppression{};
        PetscCallMPIAbort(
            backend_.communicator(),
            MPI_Allreduce(localPhaseSuppression.data(), globalPhaseSuppression.data(),
                          static_cast<int>(localPhaseSuppression.size()), MPI_LONG_LONG, MPI_SUM,
                          backend_.communicator()));

        result.bulkVolume = globalScalarSum[0];
        result.poreVolume = globalScalarSum[1];
        result.pressureAverage = result.bulkVolume > 0.0 ? globalScalarSum[2] / result.bulkVolume : 0.0;
        result.porosityAverage = result.bulkVolume > 0.0 ? globalScalarSum[3] / result.bulkVolume : 0.0;
        result.pressureMinimum = result.cellCount > 0 ? globalMin[0] : 0.0;
        result.porosityMinimum = result.cellCount > 0 ? globalMin[1] : 0.0;
        result.pressureMaximum = result.cellCount > 0 ? globalMax[0] : 0.0;
        result.porosityMaximum = result.cellCount > 0 ? globalMax[1] : 0.0;
        result.aqueousCO2MassFractionAverage =
            globalScalarSum[4] > 0.0 ? globalScalarSum[5] / globalScalarSum[4] : 0.0;
        result.aqueousCO2MassFractionMaximum = globalMax[2];
        result.trappedGasSaturationAverage =
            result.poreVolume > 0.0 ? globalScalarSum[6] / result.poreVolume : 0.0;
        result.trappedGasSaturationMaximum = globalMax[3];

        for (std::size_t phase = 0; phase < phaseCount; ++phase)
        {
            result.phasePoreVolume[phase] = globalPhaseVolume[phase];
            result.saturationAverage[phase] =
                result.poreVolume > 0.0 ? globalSaturationSum[phase] / result.poreVolume : 0.0;
            result.saturationMinimum[phase] = result.cellCount > 0 ? globalSaturationMin[phase] : 0.0;
            result.saturationMaximum[phase] = result.cellCount > 0 ? globalSaturationMax[phase] : 0.0;

            if (globalPhaseVolume[phase] > 0.0)
            {
                result.densityAverage[phase] = globalDensitySum[phase] / globalPhaseVolume[phase];
                result.viscosityAverage[phase] = globalViscositySum[phase] / globalPhaseVolume[phase];
                result.densityMinimum[phase] = globalDensityMin[phase];
                result.densityMaximum[phase] = globalDensityMax[phase];
                result.viscosityMinimum[phase] = globalViscosityMin[phase];
                result.viscosityMaximum[phase] = globalViscosityMax[phase];
            }
        }

        result.liquidOnlyCells = globalStates[0];
        result.vaporOnlyCells = globalStates[1];
        result.twoPhaseCells = globalStates[2];
        result.phasePresenceCells = globalPhasePresence;
        result.phaseSuppressionCells = globalPhaseSuppression;
        return result;
    }

