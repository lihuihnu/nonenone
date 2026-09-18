/**
 * @file natural_petsc_runtime_common.inc
 * @brief Natural/PETSc 运行时通用校验、Vec 访问和状态编解码辅助实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    /**
     * @brief Query PETSc Vec object state when the public API exists.
     *
     * PETSc 3.19 (Ubuntu 24.04 distro package used by CI) exposes neither
     * VecGetState nor PetscObjectStateGet publicly.  On those older versions
     * return a fresh synthetic state on every query.  This deliberately
     * disables cache reuse while preserving correctness; no private PETSc
     * headers or object-layout assumptions are used.
     */
    PetscErrorCode vectorState_(Vec vector, PetscObjectState *state) const noexcept
    {
#if PETSC_VERSION_GE(3,24,0)
        return VecGetState(vector, state);
#else
        (void)vector;
        *state = ++compatibilityVectorStateCounter_;
        return 0;
#endif
    }

    void validateOptions_() const
    {
        if (!(options_.timeStep > 0.0) ||
            !std::isfinite(options_.timeStep))
        {
            throw std::invalid_argument(
                "Natural time step must be finite and positive.");
        }

        if (!(options_.fugacityScalingFactor > 0.0) ||
            !std::isfinite(
                options_.fugacityScalingFactor))
        {
            throw std::invalid_argument(
                "Fugacity scaling factor must be finite and positive.");
        }

        if (!(options_.landConstant >= 0.0) ||
            !std::isfinite(
                options_.landConstant))
        {
            throw std::invalid_argument(
                "Land constant must be finite and non-negative.");
        }

        if (!(options_.rockDensity > 0.0) ||
            !std::isfinite(
                options_.rockDensity))
        {
            throw std::invalid_argument(
                "Rock density must be finite and positive.");
        }

        options_.scaling.validate();

        if (options_.maximumAqueousUnsupportedNewtonMoleFraction)
        {
            if constexpr (!Indices::fullyCompositionalThreePhase)
            {
                throw std::invalid_argument(
                    "Aqueous Newton composition limiting requires a fully compositional three-phase model.");
            }
            const double maximum =
                *options_.maximumAqueousUnsupportedNewtonMoleFraction;
            if (!(maximum >= 0.0) || !std::isfinite(maximum) ||
                options_.aqueousNewtonWaterComponent < 0 ||
                options_.aqueousNewtonWaterComponent >= Indices::numComponents ||
                options_.aqueousNewtonCo2Component < 0 ||
                options_.aqueousNewtonCo2Component >= Indices::numComponents ||
                options_.aqueousNewtonWaterComponent == options_.aqueousNewtonCo2Component)
            {
                throw std::invalid_argument(
                    "Invalid aqueous Newton composition-domain configuration.");
            }
        }

        if constexpr (Indices::hasAdsorption)
        {
            if (!options_.standardGasDensity)
            {
                throw std::invalid_argument(
                    "Adsorption runtime requires component standard-gas densities.");
            }
        }
    }

    void requireReadyForResidual_() const
    {
        requirePhaseStateInitialized_();

        if (!historyInitialized_)
        {
            throw std::logic_error(
                "Natural previous-time state is not initialized. Call initializeHistory().");
        }

        if (!backend_.isSetup())
        {
            throw std::logic_error(
                "Grid setup() must be completed before Natural residual assembly.");
        }

        if (previousPropertiesCache_.size() != localCacheBlockCount_ ||
            previousFluidAccumulationCache_.size() != localCacheBlockCount_ ||
            (Indices::hasAdsorption &&
             previousAdsorbedAccumulationCache_.size() != localCacheBlockCount_))
        {
            throw std::logic_error(
                "Natural accepted-history property/accumulation cache is incomplete.");
        }
    }

    void requirePhaseStateInitialized_() const
    {
        if (!phaseStateInitialized_)
        {
            throw std::logic_error(
                "Natural phase-state is not initialized.");
        }
    }

    static void requireVector_(
        Vec vector,
        const char *label)
    {
        if (vector == nullptr)
        {
            throw std::invalid_argument(
                std::string(label) +
                " Vec is null.");
        }
    }

    static void destroyVector_(
        Vec &vector) noexcept
    {
        if (vector != nullptr)
        {
            PetscCallAbort(
                PetscObjectComm(
                    reinterpret_cast<PetscObject>(
                        vector)),
                VecDestroy(&vector));
        }
    }

    [[nodiscard]] int resolveCO2Component_() const
    {
        if constexpr (!Indices::hasAqueousCO2Dissolution)
        {
            return -1;
        }
        else
        {
            const int component =
                options_.dissolvedCO2Component >= 0
                    ? options_.dissolvedCO2Component
                    : fluid_.aqueousCO2Component;

            if (component < 0 ||
                component >= Indices::numComponents)
            {
                throw std::invalid_argument(
                    "Invalid dissolved CO2 component index.");
            }

            return component;
        }
    }

    /** @brief 返回构造阶段已经校验并解析好的 CO2 组分号。 */
    [[nodiscard]] int resolvedCO2Component_() const noexcept
    {
        return resolvedCO2ComponentCache_;
    }

    [[nodiscard]] const std::array<
        double,
        Indices::numComponents> *
    standardGasDensityPointer_() const
    {
        if constexpr (Indices::hasAdsorption)
        {
            return &options_.standardGasDensity.value();
        }
        else
        {
            return nullptr;
        }
    }

    [[nodiscard]] static PrimaryArray
    readOwnedPrimary_(
        const PetscScalar *array,
        PetscInt begin,
        PetscInt end,
        const DofMap &map,
        const Cell &cell)
    {
        PrimaryArray result{};

        for (int component = 0;
             component <
                 Indices::numPrimaryVariables;
             ++component)
        {
            const PetscInt global =
                map.globalIndex(
                    cell,
                    component);

            if (global < begin ||
                global >= end)
            {
                throw std::runtime_error(
                    "Attempted owned primary access outside Vec ownership range.");
            }

            result[
                static_cast<std::size_t>(
                    component)] =
                PetscRealPart(
                    array[global - begin]);
        }

        return result;
    }

    static void writeOwnedPrimary_(
        PetscScalar *array,
        PetscInt begin,
        PetscInt end,
        const DofMap &map,
        const Cell &cell,
        const PrimaryArray &values)
    {
        for (int component = 0;
             component <
                 Indices::numPrimaryVariables;
             ++component)
        {
            const PetscInt global =
                map.globalIndex(
                    cell,
                    component);

            if (global < begin ||
                global >= end)
            {
                throw std::runtime_error(
                    "Attempted owned primary write outside Vec ownership range.");
            }

            array[global - begin] =
                static_cast<PetscScalar>(
                    values[
                        static_cast<std::size_t>(
                            component)]);
        }
    }

    [[nodiscard]] static PhaseStateData<Indices>
    readOwnedPhaseState_(
        const PetscScalar *array,
        PetscInt begin,
        PetscInt end,
        const DofMap &map,
        const Cell &cell)
    {
        typename PetscPhaseStateCodec<Indices>::Array
            values{};

        for (int component = 0;
             component <
                 Indices::numPhaseStateVariables;
             ++component)
        {
            const PetscInt global =
                map.globalIndex(
                    cell,
                    component);

            if (global < begin ||
                global >= end)
            {
                throw std::runtime_error(
                    "Attempted owned phase-state access outside Vec ownership range.");
            }

            values[
                static_cast<std::size_t>(
                    component)] =
                PetscRealPart(
                    array[global - begin]);
        }

        return PetscPhaseStateCodec<Indices>::decode(
            values);
    }

    static void writeOwnedPhaseState_(
        PetscScalar *array,
        PetscInt begin,
        PetscInt end,
        const DofMap &map,
        const Cell &cell,
        const PhaseStateData<Indices> &state)
    {
        const auto values =
            PetscPhaseStateCodec<Indices>::encode(
                state);

        for (int component = 0;
             component <
                 Indices::numPhaseStateVariables;
             ++component)
        {
            const PetscInt global =
                map.globalIndex(
                    cell,
                    component);

            if (global < begin ||
                global >= end)
            {
                throw std::runtime_error(
                    "Attempted owned phase-state write outside Vec ownership range.");
            }

            array[global - begin] =
                static_cast<PetscScalar>(
                    values[
                        static_cast<std::size_t>(
                            component)]);
        }
    }

    template <class ResidualBlock>
    void writeOwnedResidualBlock_(
        PetscScalar *array,
        PetscInt begin,
        PetscInt end,
        PetscInt globalBlock,
        PetscInt localBlock,
        const ResidualBlock &block)
    {
        constexpr PetscInt N =
            static_cast<PetscInt>(Indices::numEquations);
        const PetscInt firstGlobal = globalBlock * N;

        if (firstGlobal < begin ||
            firstGlobal + N > end)
        {
            throw std::runtime_error(
                "Attempted residual write outside Vec ownership range.");
        }

        for (PetscInt equation = 0; equation < N; ++equation)
        {
            const double scale = options_.scaling.enabled
                ? equationScaleFactor_(equation, localBlock)
                : 1.0;
            array[firstGlobal + equation - begin] =
                static_cast<PetscScalar>(
                    scalarValue(
                        block.value[
                            static_cast<std::size_t>(
                                equation)]) * scale);
        }
    }

    template <class LocalView>
    [[nodiscard]] static PrimaryArray
    readPrimaryLocal_(
        const LocalView &view,
        PetscInt localBlock)
    {
        PrimaryArray result{};
        for (int component = 0;
             component < Indices::numPrimaryVariables;
             ++component)
        {
            result[static_cast<std::size_t>(component)] =
                view.valueByLocalBlock(localBlock, component);
        }
        return result;
    }

    template <class LocalView>
    [[nodiscard]] static PhaseStateData<Indices>
    readPhaseStateLocal_(
        const LocalView &view,
        PetscInt localBlock)
    {
        typename PetscPhaseStateCodec<Indices>::Array values{};
        for (int component = 0;
             component < Indices::numPhaseStateVariables;
             ++component)
        {
            values[static_cast<std::size_t>(component)] =
                view.valueByLocalBlock(localBlock, component);
        }
        return PetscPhaseStateCodec<Indices>::decode(values);
    }

