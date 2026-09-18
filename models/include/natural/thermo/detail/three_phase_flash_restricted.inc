/**
 * @file three_phase_flash_restricted.inc
 * @brief 三相 Flash restricted 单/两相求解、恢复和稳定 reduced-set 后备实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    [[nodiscard]] Result flashRestrictedImpl_(
        double pressure,
        double temperature,
        Composition z,
        PhasePresence allowed,
        const std::array<Composition, 3> *phaseCompositionSeed) const
    {
        validatePT_(pressure, temperature);
        normalize_(z);
        if (allowed.empty())
            throw std::invalid_argument("Restricted flash requires at least one allowed phase.");

        if (allowed.count() == 3)
            return flash(pressure, temperature, z);
        if (allowed.count() == 1)
            return singlePhaseResult_(pressure, temperature, z, allowed);

        const auto pair = activePair_(allowed);
        Composition K = pairInitialK_(
            pair.first, pair.second, pressure, temperature, z);

        if (phaseCompositionSeed != nullptr)
        {
            Composition first = phaseCompositionSeed->at(
                static_cast<std::size_t>(phaseIndex(pair.first)));
            Composition second = phaseCompositionSeed->at(
                static_cast<std::size_t>(phaseIndex(pair.second)));
            normalize_(first);
            normalize_(second);
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                K[c] = std::clamp(
                    second[c] / std::max(first[c], options_.compositionFloor),
                    1.0e-100, 1.0e100);
            }
        }

        Result result = flashTwo_(
            pressure, temperature, z, pair.first, pair.second, K);
        if (!result.converged && phaseCompositionSeed != nullptr)
        {
            // Fall back to the generic pair seed if the current Newton
            // compositions were too far from the reduced equilibrium basin.
            K = pairInitialK_(
                pair.first, pair.second, pressure, temperature, z);
            result = flashTwo_(
                pressure, temperature, z, pair.first, pair.second, K);
        }

        if (result.converged)
            canonicalize_(pressure, temperature, result);
        return result;
    }

    struct TwoPhaseLogKState
    {
        bool valid{false};
        double norm{std::numeric_limits<double>::infinity()};
        Composition residual{};
        Result result{};
    };

    [[nodiscard]] TwoPhaseLogKState evaluateTwoPhaseLogK_(
        double pressure,
        double temperature,
        const Composition &z,
        CompositionalPhase first,
        CompositionalPhase second,
        const Composition &logK,
        double preferredSecond) const
    {
        TwoPhaseLogKState state;
        Composition K{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            K[c] = std::exp(std::clamp(logK[c], -230.0, 230.0));
        }
        const Rr2 rr = solveRr2_(z, K, preferredSecond);
        if (!rr.converged)
            return state;
        const double betaSecond = rr.betaSecond;
        const double betaFirst = 1.0 - betaSecond;

        Composition a{}, b{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double D = (1.0 - betaSecond) + betaSecond * K[c];
            if (!(D > 0.0) || !std::isfinite(D))
                return state;
            a[c] = z[c] / D;
            b[c] = K[c] * a[c];
        }
        normalize_(a);
        normalize_(b);
        const auto ta = phaseResult_(first, pressure, temperature, a);
        const auto tb = phaseResult_(second, pressure, temperature, b);

        state.norm = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double phiA = std::max(ta.fugacityCoefficient[c], 1.0e-300);
            const double phiB = std::max(tb.fugacityCoefficient[c], 1.0e-300);
            state.residual[c] = logK[c] - std::log(phiA / phiB);
            state.norm = std::max(state.norm, std::abs(state.residual[c]));
        }

        state.result.presence = PhasePresence(
            PhasePresence::bit(first) | PhasePresence::bit(second));
        state.result.phaseMoleFraction.fill(0.0);
        state.result.phaseMoleFraction[static_cast<std::size_t>(phaseIndex(first))] = betaFirst;
        state.result.phaseMoleFraction[static_cast<std::size_t>(phaseIndex(second))] = betaSecond;
        state.result.composition[static_cast<std::size_t>(phaseIndex(first))] = a;
        state.result.composition[static_cast<std::size_t>(phaseIndex(second))] = b;
        state.result.compressibility[static_cast<std::size_t>(phaseIndex(first))] = ta.compressibility;
        state.result.compressibility[static_cast<std::size_t>(phaseIndex(second))] = tb.compressibility;
        for (CompositionalPhase phase : phases_)
        {
            if (!state.result.presence.contains(phase))
                state.result.composition[static_cast<std::size_t>(phaseIndex(phase))] = z;
        }
        state.valid = std::isfinite(state.norm);
        return state;
    }

    [[nodiscard]] Result flashTwo_(
        double pressure,
        double temperature,
        const Composition &z,
        CompositionalPhase first,
        CompositionalPhase second,
        Composition K) const
    {
        Result result;
        result.presence = PhasePresence(
            PhasePresence::bit(first) | PhasePresence::bit(second));
        for (double &value : K)
            value = std::clamp(value, 1.0e-100, 1.0e100);

        double betaSecond = 0.5;
        int totalIterations = 0;
        const int ssiSteps = std::min(
            options_.maximumIterations,
            std::max(options_.twoPhaseSsiIterationsBeforeNewton, 0));

        // 数值：先用阻尼 SSI 把 K 值推入正确吸引域，再切换到 log-K Newton 加速局部收敛。
        auto oneSsiStep = [&]() -> bool {
            const Rr2 rr = solveRr2_(z, K, betaSecond);
            if (!rr.converged)
                return false;
            betaSecond = rr.betaSecond;
            const double betaFirst = 1.0 - betaSecond;
            Composition a{}, b{};
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double D = (1.0 - betaSecond) + betaSecond * K[c];
                if (!(D > 0.0) || !std::isfinite(D))
                    return false;
                a[c] = z[c] / D;
                b[c] = K[c] * a[c];
            }
            normalize_(a);
            normalize_(b);
            const auto ta = phaseResult_(first, pressure, temperature, a);
            const auto tb = phaseResult_(second, pressure, temperature, b);

            double maxResidual = 0.0;
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double target =
                    std::max(ta.fugacityCoefficient[c], 1.0e-300) /
                    std::max(tb.fugacityCoefficient[c], 1.0e-300);
                const double d = std::clamp(
                    std::log(target / K[c]),
                    -options_.logKStepLimit, options_.logKStepLimit);
                maxResidual = std::max(maxResidual, std::abs(d));
                K[c] *= std::exp(0.6 * d);
                K[c] = std::clamp(K[c], 1.0e-100, 1.0e100);
            }

            ++totalIterations;
            result.iterations = totalIterations;
            result.phaseMoleFraction.fill(0.0);
            result.phaseMoleFraction[static_cast<std::size_t>(phaseIndex(first))] = betaFirst;
            result.phaseMoleFraction[static_cast<std::size_t>(phaseIndex(second))] = betaSecond;
            result.composition[static_cast<std::size_t>(phaseIndex(first))] = a;
            result.composition[static_cast<std::size_t>(phaseIndex(second))] = b;
            result.compressibility[static_cast<std::size_t>(phaseIndex(first))] = ta.compressibility;
            result.compressibility[static_cast<std::size_t>(phaseIndex(second))] = tb.compressibility;
            for (CompositionalPhase phase : phases_)
            {
                if (!result.presence.contains(phase))
                    result.composition[static_cast<std::size_t>(phaseIndex(phase))] = z;
            }

            if (maxResidual < options_.fugacityTolerance)
            {
                result.converged = true;
                result.presence = presenceFromFractions_(result.phaseMoleFraction);
                deriveReferenceK_(result);
                updateSaturationsFromMoles(pressure, temperature, result);
            }
            return true;
        };

        for (int iteration = 0; iteration < ssiSteps; ++iteration)
        {
            if (!oneSsiStep())
                return result;
            if (result.converged)
                return result;
        }

        Composition logK{};
        for (int i = 0; i < N; ++i)
            logK[static_cast<std::size_t>(i)] =
                std::log(K[static_cast<std::size_t>(i)]);

        for (int iteration = 0;
             iteration < options_.maximumTwoPhaseNewtonIterations;
             ++iteration)
        {
            const auto base = evaluateTwoPhaseLogK_(
                pressure, temperature, z, first, second, logK, betaSecond);
            if (!base.valid)
                break;
            result = base.result;
            result.iterations = ++totalIterations;
            betaSecond = result.phaseMoleFraction[static_cast<std::size_t>(phaseIndex(second))];
            if (base.norm < options_.fugacityTolerance)
            {
                result.converged = true;
                result.presence = presenceFromFractions_(result.phaseMoleFraction);
                deriveReferenceK_(result);
                updateSaturationsFromMoles(pressure, temperature, result);
                return result;
            }

            std::array<std::array<double, N>, N> jacobian{};
            bool jacobianValid = true;
            for (int column = 0; column < N; ++column)
            {
                const std::size_t j = static_cast<std::size_t>(column);
                Composition perturbed = logK;
                const double step = options_.numericalJacobianStep *
                    std::max(1.0, std::abs(logK[j]));
                perturbed[j] += step;
                const auto shifted = evaluateTwoPhaseLogK_(
                    pressure, temperature, z, first, second,
                    perturbed, betaSecond);
                if (!shifted.valid)
                {
                    jacobianValid = false;
                    break;
                }
                for (int row = 0; row < N; ++row)
                {
                    const std::size_t i = static_cast<std::size_t>(row);
                    jacobian[i][j] =
                        (shifted.residual[i] - base.residual[i]) / step;
                }
            }
            if (!jacobianValid)
                break;

            Composition rhs{}, delta{};
            for (int i = 0; i < N; ++i)
                rhs[static_cast<std::size_t>(i)] =
                    -base.residual[static_cast<std::size_t>(i)];
            if (!solveDenseLinearSystem_<static_cast<std::size_t>(N)>(
                    jacobian, rhs, delta))
                break;

            bool accepted = false;
            double damping = 1.0;
            for (int lineSearch = 0; lineSearch < 16; ++lineSearch)
            {
                Composition trialLogK = logK;
                for (int i = 0; i < N; ++i)
                {
                    const std::size_t c = static_cast<std::size_t>(i);
                    trialLogK[c] += damping * delta[c];
                }
                const auto trial = evaluateTwoPhaseLogK_(
                    pressure, temperature, z, first, second,
                    trialLogK, betaSecond);
                if (trial.valid && trial.norm < base.norm)
                {
                    logK = trialLogK;
                    accepted = true;
                    break;
                }
                damping *= 0.5;
            }
            if (!accepted)
                break;
        }

        for (int iteration = ssiSteps;
             iteration < options_.maximumIterations;
             ++iteration)
        {
            for (int i = 0; i < N; ++i)
                K[static_cast<std::size_t>(i)] = std::exp(std::clamp(
                    logK[static_cast<std::size_t>(i)], -230.0, 230.0));
            if (!oneSsiStep())
                return result;
            if (result.converged)
                return result;
            for (int i = 0; i < N; ++i)
                logK[static_cast<std::size_t>(i)] =
                    std::log(K[static_cast<std::size_t>(i)]);
        }
        return result;
    }

    /**
     * @brief 终检一个程序路径标记失败、但物理残差已满足容差的两相结果。
     *
     * 仅在 restricted flash 的常规求解已返回 `converged=false` 后调用。
     * 接受条件同时要求：两相分率有效且闭合、总体组成守恒、逐组分逸度相等。
     * 因而这不是放宽 flash 方程，而是用真正的物理方程校正内部流程标志。
     */
    void finalizePhysicallyConvergedTwoPhase_(
        double pressure,
        double temperature,
        const Composition &z,
        CompositionalPhase first,
        CompositionalPhase second,
        Result &result) const
    {
        if (result.converged)
            return;

        const std::size_t a = static_cast<std::size_t>(phaseIndex(first));
        const std::size_t b = static_cast<std::size_t>(phaseIndex(second));
        const double betaA = result.phaseMoleFraction[a];
        const double betaB = result.phaseMoleFraction[b];

        // Do not reinterpret a collapsed single-phase endpoint as a two-phase
        // equilibrium; existing single-phase active-set logic handles it.
        if (!(betaA > options_.phaseFractionTolerance) ||
            !(betaB > options_.phaseFractionTolerance) ||
            !std::isfinite(betaA) || !std::isfinite(betaB))
            return;

        const double betaClosure = std::abs(betaA + betaB - 1.0);
        const double closureTolerance = std::max(
            options_.fugacityTolerance,
            10.0 * options_.rachfordRiceTolerance);
        if (betaClosure > closureTolerance)
            return;

        Composition xa = result.composition[a];
        Composition xb = result.composition[b];
        normalize_(xa);
        normalize_(xb);

        double materialClosure = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double reconstructed = betaA * xa[c] + betaB * xb[c];
            materialClosure = std::max(
                materialClosure,
                std::abs(reconstructed - z[c]));
        }
        if (materialClosure > closureTolerance)
            return;

        const auto thermoA = phaseResult_(
            first, pressure, temperature, xa);
        const auto thermoB = phaseResult_(
            second, pressure, temperature, xb);

        double fugacityResidual = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double activityA = std::max(
                xa[c] * thermoA.fugacityCoefficient[c], 1.0e-300);
            const double activityB = std::max(
                xb[c] * thermoB.fugacityCoefficient[c], 1.0e-300);
            fugacityResidual = std::max(
                fugacityResidual,
                std::abs(std::log(activityA / activityB)));
        }
        if (fugacityResidual > options_.fugacityTolerance)
            return;

        result.composition[a] = xa;
        result.composition[b] = xb;
        result.presence = PhasePresence(
            PhasePresence::bit(first) | PhasePresence::bit(second));
        result.compressibility[a] = thermoA.compressibility;
        result.compressibility[b] = thermoB.compressibility;
        result.converged = true;
        deriveReferenceK_(result);
        updateSaturationsFromMoles(pressure, temperature, result);
    }

    void deriveReferenceK_(Result &result) const
    {
        const Composition &xo = result.composition[0];
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            result.vaporOilK[c] = result.composition[1][c] / std::max(xo[c], options_.compositionFloor);
            result.waterOilK[c] = result.composition[2][c] / std::max(xo[c], options_.compositionFloor);
        }
    }

    [[nodiscard]] Result singlePhaseResult_(
        double pressure,
        double temperature,
        const Composition &z,
        PhasePresence active) const
    {
        Result result;
        result.converged = true;
        result.presence = active;
        result.phaseMoleFraction.fill(0.0);
        const CompositionalPhase p = firstActive_(active);
        const int pi = phaseIndex(p);
        result.phaseMoleFraction[static_cast<std::size_t>(pi)] = 1.0;
        result.composition[static_cast<std::size_t>(pi)] = z;
        for (CompositionalPhase q : phases_)
            if (!active.contains(q))
                result.composition[static_cast<std::size_t>(phaseIndex(q))] = z;
        const auto thermo = phaseResult_(p, pressure, temperature, z);
        result.compressibility[static_cast<std::size_t>(pi)] = thermo.compressibility;
        deriveReferenceK_(result);
        updateSaturationsFromMoles(pressure, temperature, result);
        return result;
    }

    [[nodiscard]] PhasePresence chooseSinglePhase_(
        double pressure,
        double temperature,
        const Composition &z) const
    {
        // 同一总体组成 z 下比较三个固定热力学相身份的 Gibbs。PR 中 O/W
        // 通常退化为同一个 liquid-root 模型，此时仅用 water-dominant 作为标签 tie-break；
        // SW 中 aqueous BIP 与 non-aqueous BIP 不同，因此 W 必须真正参加 Gibbs 比较。
        const auto oil = phaseResult_(CompositionalPhase::Oil, pressure, temperature, z);
        const auto gas = phaseResult_(CompositionalPhase::Gas, pressure, temperature, z);
        const auto water = phaseResult_(CompositionalPhase::Water, pressure, temperature, z);
        double go = 0.0, gg = 0.0, gw = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            go += z[c] * std::log(std::max(oil.fugacityCoefficient[c], 1.0e-300));
            gg += z[c] * std::log(std::max(gas.fugacityCoefficient[c], 1.0e-300));
            gw += z[c] * std::log(std::max(water.fugacityCoefficient[c], 1.0e-300));
        }
        constexpr double tieTolerance = 1.0e-10;
        if (gg + tieTolerance < go && gg + tieTolerance < gw)
            return PhasePresence::gasOnly();
        if (gw + tieTolerance < go && gw + tieTolerance < gg)
            return PhasePresence::waterOnly();
        if (waterIsDominant_(z) && std::abs(gw - go) <= tieTolerance)
            return PhasePresence::waterOnly();
        if (gg + tieTolerance < go)
            return PhasePresence::gasOnly();
        return PhasePresence::oilOnly();
    }

    [[nodiscard]] double dimensionlessGibbs_(
        double pressure,
        double temperature,
        const Result &result) const
    {
        double g = 0.0;
        for (CompositionalPhase phase : phases_)
        {
            const std::size_t slot = static_cast<std::size_t>(phaseIndex(phase));
            const double beta = result.phaseMoleFraction[slot];
            if (!(beta > options_.phaseFractionTolerance))
                continue;
            const auto thermo = phaseResultForThermodynamicRole_(
                phase, pressure, temperature, result.composition[slot]);
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double x = std::max(
                    result.composition[slot][c], options_.compositionFloor);
                const double phi = std::max(
                    thermo.fugacityCoefficient[c], 1.0e-300);
                g += beta * result.composition[slot][c] *
                    (std::log(x) + std::log(phi));
            }
        }
        return g;
    }

    /**
     * @brief 将公开标准 O/G/W 结果逆映射回当前 SW 求解器的固定热力学角色槽位。
     *
     * `canonicalizeSwThermodynamicRoles_()` 在对外边界执行“数值槽位 -> 物理 O/G/W”映射；
     * 压力延拓作为内部初值时必须反向映射，否则会把 Aqueous BIP 或气/液根附到错误组成。
     */
    [[nodiscard]] Result numericalSwRoleResult_(const Result &canonical) const
    {
        if (!eos_.usesSoreideWhitson())
            return canonical;

        Result numerical = canonical;
        numerical.phaseMoleFraction = {0.0, 0.0, 0.0};
        numerical.saturation = {0.0, 0.0, 0.0};
        numerical.compressibility = {1.0, 1.0, 1.0};
        numerical.molarDensity = {0.0, 0.0, 0.0};
        numerical.composition = {};

        for (int numericalSlot = 0; numericalSlot < 3; ++numericalSlot)
        {
            const std::size_t destination =
                static_cast<std::size_t>(numericalSlot);
            const std::size_t physical = static_cast<std::size_t>(
                phaseIndex(thermodynamicRoles_[destination]));
            numerical.phaseMoleFraction[destination] =
                canonical.phaseMoleFraction[physical];
            numerical.saturation[destination] = canonical.saturation[physical];
            numerical.composition[destination] = canonical.composition[physical];
            numerical.compressibility[destination] =
                canonical.compressibility[physical];
            numerical.molarDensity[destination] = canonical.molarDensity[physical];
        }
        numerical.presence = presenceFromFractions_(numerical.phaseMoleFraction);
        deriveReferenceK_(numerical);
        return numerical;
    }

    /**
     * @brief 仅在所有常规路径失败后启用的通用压力同伦恢复。
     *
     * 数值：PR/SW/CPA 都可能在低温、临界附近或窄三相区出现 Wilson 初值与真实平衡支路
     * 之间的吸引域间隙，而邻近压力仍可稳定收敛。算法先寻找同一 EOS/固定热力学角色下
     * 可收敛的邻近压力锚点，再沿对数压力小步延拓到目标压力。每一步仍调用正常
     * active-set/稳定性逻辑，因此相出现/消失继续由热力学判据控制；它只改变失败恢复初值，
     * 不改变任何已收敛状态的热力学判据。
     */
    [[nodiscard]] Result pressureContinuationFallback_(
        double pressure,
        double temperature,
        const Composition &z) const
    {
        Result failure;
        if (!(pressure > 0.0))
            return failure;

        constexpr std::array<double, 8> anchorRatios{{
            0.95, 1.05, 0.90, 1.10, 0.80, 1.25, 0.65, 1.50}};
        // The Jia boundary audit can contain very narrow CPA three-phase
        // basins. Once the ordinary flash has already failed, use a finer
        // fail-only pressure homotopy for the opt-in CPA candidate-Gibbs
        // policy. Production/default CPA, PR and SW retain the historical
        // continuation spacing and iteration cap.
        const bool fineCpaBoundaryContinuation =
            eos_.usesCubicPlusAssociation() &&
            options_.cpaSelectGibbsMinimumRoot;
        const double maximumLogPressureStep =
            fineCpaBoundaryContinuation ? 0.005 : 0.025;
        const int maximumContinuationSteps =
            fineCpaBoundaryContinuation ? 40 : 12;

        // 数值：锚点可使用既有恢复分支，但禁止再次启动压力延拓或角色枚举，避免递归搜索。
        CubicThreePhaseFlash anchorSolver(
            eos_, options_, thermodynamicRoles_, false, false);

        const double logTarget = std::log(pressure);
        for (const double ratio : anchorRatios)
        {
            const double anchorPressure = pressure * ratio;
            if (!(anchorPressure > 0.0) || !std::isfinite(anchorPressure))
                continue;

            Result anchor = anchorSolver.flash(
                anchorPressure, temperature, z);
            if (!anchor.converged)
                continue;

            Result seed = eos_.usesSoreideWhitson()
                ? numericalSwRoleResult_(anchor)
                : anchor;
            if (!seed.converged || seed.presence.empty())
                continue;

            const double logAnchor = std::log(anchorPressure);
            const double logSpan = std::abs(logTarget - logAnchor);
            const int continuationSteps = std::clamp(
                static_cast<int>(std::ceil(logSpan / maximumLogPressureStep)),
                1, maximumContinuationSteps);
            int totalIterations = anchor.iterations;
            bool pathConverged = true;
            for (int step = 1; step <= continuationSteps; ++step)
            {
                const double fraction =
                    static_cast<double>(step) / continuationSteps;
                const double stepPressure = std::exp(
                    logAnchor + fraction * (logTarget - logAnchor));

                Result current = equilibrateActiveSet_(
                    stepPressure, temperature, z, seed.presence, seed, false);
                if (!current.converged)
                {
                    pathConverged = false;
                    break;
                }
                totalIterations += current.iterations;
                seed = current;
            }

            if (pathConverged && seed.converged)
            {
                seed.iterations = totalIterations;
                return seed;
            }
        }
        return failure;
    }

    /**
     * @brief fail-only reduced-set recovery 使用的 restricted flash。
     *
     * PR/CPA 直接复用公共 O/G/W 规范化；SW 外部 active-set 明确指定 G+W/O+W 时
     * 可在 cubic 根合并区保留该角色，但 recovery 枚举多个候选相集时仍通过历史
     * 公共角色检查，避免改变已有 CO2-brine 和 P-T 扫描的候选排序。
     */
    [[nodiscard]] Result flashRestrictedRecovery_(
        double pressure,
        double temperature,
        Composition z,
        PhasePresence allowed) const
    {
        normalize_(z);
        Result result = flashRestrictedImpl_(
            pressure, temperature, z, allowed, nullptr);
        if (result.converged)
        {
            canonicalizeFinal_(pressure, temperature, result, false);
            if (result.converged &&
                !validatePhysicalFlashResult_(pressure, temperature, z, result))
                result.converged = false;
        }
        return result;
    }

    /**
     * @brief 通用三 EOS 失败恢复：搜索稳定的一/两相约化解。
     *
     * 数值：相边界附近，TPD 可能提示第三相出现，但直接三相 Newton/SSI 也可能
     * 没有进入正确吸引域。枚举 O+G、O+W、G+W，候选必须通过完整缺失相
     * stability.valid && stability.stable 证书，并按 Gibbs 最低者选择。
     */
    [[nodiscard]] Result stableReducedSetFallback_(
        double pressure,
        double temperature,
        const Composition &z) const
    {
        Result best;
        double bestG = std::numeric_limits<double>::infinity();
        const std::array<PhasePresence, 3> pairs{{
            PhasePresence(PhasePresence::bit(CompositionalPhase::Oil) |
                          PhasePresence::bit(CompositionalPhase::Gas)),
            PhasePresence(PhasePresence::bit(CompositionalPhase::Oil) |
                          PhasePresence::bit(CompositionalPhase::Water)),
            PhasePresence(PhasePresence::bit(CompositionalPhase::Gas) |
                          PhasePresence::bit(CompositionalPhase::Water))}};

        for (const PhasePresence allowed : pairs)
        {
            Result candidate = flashRestrictedRecovery_(
                pressure, temperature, z, allowed);
            if (!candidate.converged)
                continue;
            candidate.presence = presenceFromFractions_(
                candidate.phaseMoleFraction);
            if (candidate.presence.empty())
                continue;
            const StabilityResult stability = stabilityTest(
                pressure, temperature, z,
                candidate.presence, candidate.composition);
            if (!stability.valid)
                continue;
            if (!stability.stable)
            {
                // A converged liquid branch is a seed, not an equilibrium
                // certificate. Release its unstable missing phase and solve
                // again with fixed active-root identities before comparing G.
                if (!eos_.usesCubicPlusAssociation() ||
                    !options_.cpaSelectGibbsMinimumRoot)
                    continue;
                candidate = equilibrateActiveSet_(
                    pressure, temperature, z, candidate.presence, candidate);
                if (!candidate.converged ||
                    !validatePhysicalFlashResult_(pressure, temperature, z, candidate))
                    continue;
                const auto expandedStability = stabilityTest(
                    pressure, temperature, z, candidate.presence, candidate.composition);
                if (!expandedStability.valid || !expandedStability.stable)
                    continue;
            }
            const double g = dimensionlessGibbs_(
                pressure, temperature, candidate);
            if (std::isfinite(g) && g < bestG)
            {
                bestG = g;
                best = candidate;
            }
        }
        return best;
    }

    /**
     * @brief SW 单相候选返回前，用稳定约化相集和 Gibbs 能排除漏相局部盆地。
     *
     * SW 的水相温度修正会改变单相 TPD 目标的形状；单一 Wilson 初值偶尔会
     * 漏掉多相分裂盆地。只对已经通过物理终检的 SW 单相候选执行本复核；
     * 稳定约化解仅在 Gibbs 能严格更低时替换。若单相 TPD 明确不稳定但尚未
     * 找到稳定约化解，则令上层继续三相/角色/延拓恢复，真实单相区保持不变。
     */
    [[nodiscard]] Result preferStableSwReducedOverSingle_(
        double pressure,
        double temperature,
        const Composition &z,
        Result candidate) const
    {
        if (!eos_.usesSoreideWhitson() ||
            presenceFromFractions_(candidate.phaseMoleFraction).count() != 1)
        {
            return candidate;
        }

        candidate.presence = presenceFromFractions_(
            candidate.phaseMoleFraction);
        const StabilityResult singleStability = stabilityTest(
            pressure, temperature, z,
            candidate.presence, candidate.composition);
        const bool singleKnownUnstable =
            singleStability.valid && !singleStability.stable;

        Result reduced;
        double reducedG = std::numeric_limits<double>::infinity();
        const std::array<PhasePresence, 3> pairs{{
            PhasePresence(PhasePresence::bit(CompositionalPhase::Oil) |
                          PhasePresence::bit(CompositionalPhase::Gas)),
            PhasePresence(PhasePresence::bit(CompositionalPhase::Oil) |
                          PhasePresence::bit(CompositionalPhase::Water)),
            PhasePresence(PhasePresence::bit(CompositionalPhase::Gas) |
                          PhasePresence::bit(CompositionalPhase::Water))}};
        for (const PhasePresence allowed : pairs)
        {
            Result trial = flashRestricted(
                pressure, temperature, z, allowed);
            if (!trial.converged)
                continue;
            trial.presence = presenceFromFractions_(
                trial.phaseMoleFraction);
            if (trial.presence.empty())
                continue;
            const StabilityResult stability = stabilityTest(
                pressure, temperature, z,
                trial.presence, trial.composition);
            if (!stability.valid || !stability.stable)
                continue;
            const double trialG = dimensionlessGibbs_(
                pressure, temperature, trial);
            if (std::isfinite(trialG) && trialG < reducedG)
            {
                reducedG = trialG;
                reduced = trial;
            }
        }
        if (!reduced.converged)
        {
            if (singleKnownUnstable)
                candidate.converged = false;
            return candidate;
        }

        const double candidateG = dimensionlessGibbs_(
            pressure, temperature, candidate);
        const double gibbsScale = std::max({
            1.0, std::abs(candidateG), std::abs(reducedG)});
        constexpr double gibbsComparisonTolerance = 1.0e-12;
        if (std::isfinite(reducedG) &&
            (!std::isfinite(candidateG) ||
             reducedG < candidateG - gibbsComparisonTolerance * gibbsScale))
        {
            return reduced;
        }
        if (singleKnownUnstable)
            candidate.converged = false;
        return candidate;
    }


