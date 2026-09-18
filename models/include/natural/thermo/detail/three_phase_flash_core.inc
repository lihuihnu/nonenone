/**
 * @file three_phase_flash_core.inc
 * @brief 三相 Flash 核心相角色映射、EOS 相评价、输入规范化与初猜实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    struct Trial
    {
        bool valid{false};
        bool converged{false};
        bool trivial{false};
        double sum{0.0};
        Composition composition{};
    };

    struct Rr3
    {
        bool converged{false};
        double betaGas{0.0};
        double betaWater{0.0};
    };

    struct Rr2
    {
        bool converged{false};
        double betaSecond{0.0};
    };

    inline static constexpr std::array<CompositionalPhase, 3> phases_{
        CompositionalPhase::Oil,
        CompositionalPhase::Gas,
        CompositionalPhase::Water};

    inline static constexpr std::array<std::array<CompositionalPhase, 3>, 6>
    swRolePermutations_{{
        {{CompositionalPhase::Oil, CompositionalPhase::Gas, CompositionalPhase::Water}},
        {{CompositionalPhase::Oil, CompositionalPhase::Water, CompositionalPhase::Gas}},
        {{CompositionalPhase::Gas, CompositionalPhase::Oil, CompositionalPhase::Water}},
        {{CompositionalPhase::Gas, CompositionalPhase::Water, CompositionalPhase::Oil}},
        {{CompositionalPhase::Water, CompositionalPhase::Oil, CompositionalPhase::Gas}},
        {{CompositionalPhase::Water, CompositionalPhase::Gas, CompositionalPhase::Oil}}
    }};

    CubicThreePhaseFlash(
        const Eos &eos,
        ThreePhaseFlashOptions options,
        std::array<CompositionalPhase, 3> thermodynamicRoles,
        bool enableSwRolePermutationSearch,
        bool enablePressureContinuation)
        : CubicThreePhaseFlash(eos, std::move(options))
    {
        thermodynamicRoles_ = thermodynamicRoles;
        enableSwRolePermutationSearch_ = enableSwRolePermutationSearch;
        enablePressureContinuation_ = enablePressureContinuation;
    }

    /**
     * @brief 将公开 O/G/W 数值槽位映射到本次求解固定的热力学角色。
     *
     * SW 中 Aqueous 角色决定水相 BIP，因此该映射在一次非线性求解期间不得变化；
     * PR/CPA 不需要额外重映射，直接使用公共槽位。
     */
    [[nodiscard]] CompositionalPhase thermodynamicRole_(
        CompositionalPhase numericalSlot) const noexcept
    {
        if (!eos_.usesSoreideWhitson())
            return numericalSlot;
        return thermodynamicRoles_[static_cast<std::size_t>(phaseIndex(numericalSlot))];
    }

    struct CubicParameterCache final
    {
        bool valid{false};
        double pressure{0.0};
        double temperature{0.0};
        typename Eos::template MixingParameters<double> nonAqueous{};
        typename Eos::template MixingParameters<double> aqueous{};
    };

    void ensureCubicParameterCache_(double pressure, double temperature) const
    {
        if (cubicParameterCache_.valid &&
            cubicParameterCache_.pressure == pressure &&
            cubicParameterCache_.temperature == temperature)
        {
            return;
        }

        cubicParameterCache_.nonAqueous = eos_.mixingParameters(
            pressure, temperature, CompositionalPhase::Oil);
        cubicParameterCache_.aqueous = eos_.usesSoreideWhitson()
            ? eos_.mixingParameters(
                  pressure, temperature, CompositionalPhase::Water)
            : cubicParameterCache_.nonAqueous;
        cubicParameterCache_.pressure = pressure;
        cubicParameterCache_.temperature = temperature;
        cubicParameterCache_.valid = true;
    }

    [[nodiscard]] typename Eos::template PhaseResult<double>
    phaseResultForThermodynamicRole_(
        CompositionalPhase thermodynamicPhase,
        double pressure,
        double temperature,
        const Composition &composition) const
    {
        if (eos_.usesCubicPlusAssociation())
        {
            return eos_.phaseResult(
                pressure,
                temperature,
                composition,
                thermodynamicPhase,
                options_.cpaSelectGibbsMinimumRoot);
        }

        // 性能：一次 P–T flash 会重复计算逸度。PR 只有一套参数矩阵，SW 额外有一套水相矩阵；
        // 这里仅缓存固定 (p,T) 下不随组成变化的参数，组成混合、根选择和逸度公式保持原样。
        ensureCubicParameterCache_(pressure, temperature);
        const auto &parameters = thermodynamicPhase == CompositionalPhase::Water
            ? cubicParameterCache_.aqueous
            : cubicParameterCache_.nonAqueous;
        const auto mix = eos_.phaseMixing(composition, parameters);
        const double z = thermodynamicPhase == CompositionalPhase::Gas
            ? eos_.vaporRoot(mix.A, mix.B)
            : eos_.liquidRoot(mix.A, mix.B);
        return eos_.fugacityAtRoot(
            pressure, composition, parameters, mix, z);
    }

    [[nodiscard]] typename Eos::template PhaseResult<double>
    phaseResult_(
        CompositionalPhase phase,
        double pressure,
        double temperature,
        const Composition &composition) const
    {
        return phaseResultForThermodynamicRole_(
            thermodynamicRole_(phase), pressure, temperature, composition);
    }

    void validatePT_(double pressure, double temperature) const
    {
        if (!(pressure > 0.0) || !std::isfinite(pressure))
            throw std::invalid_argument("Three-phase flash pressure must be finite and positive.");
        if (!(temperature > 0.0) || !std::isfinite(temperature))
            throw std::invalid_argument("Three-phase flash temperature must be finite and positive.");
    }

    void normalize_(Composition &x) const
    {
        double sum = 0.0;
        for (double &v : x)
        {
            if (!std::isfinite(v))
                throw std::invalid_argument("Three-phase flash composition contains a non-finite value.");
            v = std::max(v, options_.compositionFloor);
            sum += v;
        }
        if (!(sum > 0.0))
            throw std::invalid_argument("Three-phase flash composition sum must be positive.");
        for (double &v : x) v /= sum;
    }

    /** @brief Wilson 关联式提供气/油 K 值初猜；它只用于全局化，不参与最终平衡判据。 */
    [[nodiscard]] Composition wilsonGasOilK_(double pressure, double temperature) const
    {
        Composition K{};
        for (int i = 0; i < N; ++i)
        {
            const double Tc = eos_.mixture().criticalTemperature(i);
            const double Pc = eos_.mixture().criticalPressure(i);
            const double omega = eos_.mixture().acentricFactor(i);
            const double lnK = std::log(Pc / pressure) +
                5.373 * (1.0 + omega) * (1.0 - Tc / temperature);
            K[static_cast<std::size_t>(i)] = std::clamp(std::exp(lnK), 1.0e-30, 1.0e30);
        }
        return K;
    }

    /** @brief 构造富水相/油相 K 值的通用初猜；SW 失败恢复会用更可靠的邻近状态替代它。 */
    [[nodiscard]] Composition initialWaterOilK_(const Composition &) const
    {
        Composition K{};
        K.fill(0.05);
        if (options_.waterComponent >= 0)
            K[static_cast<std::size_t>(options_.waterComponent)] = 50.0;
        else
            K.back() = 50.0;
        return K;
    }

