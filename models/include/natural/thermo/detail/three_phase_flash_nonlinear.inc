/**
 * @file three_phase_flash_nonlinear.inc
 * @brief 三相 Flash Rachford–Rice、log-K 三相非线性求解与主三相迭代实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    [[nodiscard]] Rr3 solveRr3_(
        const Composition &z,
        const Composition &kg,
        const Composition &kw,
        double initialBg,
        double initialBw) const
    {
        constexpr double eps = 1.0e-12;
        double bg = std::clamp(initialBg, 1.0e-4, 0.8);
        double bw = std::clamp(initialBw, 1.0e-4, 0.8);
        if (bg + bw >= 0.95)
        {
            const double scale = 0.8 / (bg + bw);
            bg *= scale;
            bw *= scale;
        }

        // 数学：未知量为 beta_g、beta_w；两个残差分别是广义 Rachford–Rice 方程，
        // Jacobian 在这里解析计算，外层再用线搜索保证相分率留在物理单纯形内。
        auto residualNorm = [&](double cg, double cw, double *fgOut = nullptr,
                                double *fwOut = nullptr, double *jggOut = nullptr,
                                double *jwwOut = nullptr, double *jgwOut = nullptr) {
            if (!(cg > 0.0) || !(cw > 0.0) || !(cg + cw < 1.0))
                return std::numeric_limits<double>::infinity();
            double fg = 0.0, fw = 0.0;
            double jgg = 0.0, jww = 0.0, jgw = 0.0;
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double ag = kg[c] - 1.0;
                const double aw = kw[c] - 1.0;
                const double co = 1.0 - cg - cw;
                const double D = co + cg * kg[c] + cw * kw[c];
                if (!(D > eps) || !std::isfinite(D))
                    return std::numeric_limits<double>::infinity();
                fg += z[c] * ag / D;
                fw += z[c] * aw / D;
                const double invD2 = z[c] / (D * D);
                jgg -= invD2 * ag * ag;
                jww -= invD2 * aw * aw;
                jgw -= invD2 * ag * aw;
            }
            if (fgOut) *fgOut = fg;
            if (fwOut) *fwOut = fw;
            if (jggOut) *jggOut = jgg;
            if (jwwOut) *jwwOut = jww;
            if (jgwOut) *jgwOut = jgw;
            return std::max(std::abs(fg), std::abs(fw));
        };

        for (int iteration = 0; iteration < options_.maximumRachfordRiceIterations; ++iteration)
        {
            double fg = 0.0, fw = 0.0, jgg = 0.0, jww = 0.0, jgw = 0.0;
            const double norm = residualNorm(bg, bw, &fg, &fw, &jgg, &jww, &jgw);
            if (!std::isfinite(norm))
                return {};
            if (norm < options_.rachfordRiceTolerance)
                return {true, bg, bw};

            const double det = jgg * jww - jgw * jgw;
            if (std::abs(det) < 1.0e-24 || !std::isfinite(det))
                return {};
            const double dbg = (-fg * jww + fw * jgw) / det;
            const double dbw = (-fw * jgg + fg * jgw) / det;

            double lambda = 1.0;
            bool accepted = false;
            for (int ls = 0; ls < 40; ++ls)
            {
                const double cg = bg + lambda * dbg;
                const double cw = bw + lambda * dbw;
                const double candidateNorm = residualNorm(cg, cw);
                if (std::isfinite(candidateNorm) && candidateNorm < norm)
                {
                    bg = cg;
                    bw = cw;
                    accepted = true;
                    break;
                }
                lambda *= 0.5;
            }
            if (!accepted)
                return {};
        }
        return {};
    }

    [[nodiscard]] Rr3 solveRr3Robust_(
        const Composition &z,
        const Composition &kg,
        const Composition &kw,
        double preferredGas,
        double preferredWater) const
    {
        const std::array<std::array<double, 2>, 17> guesses{{
            {{preferredGas, preferredWater}},
            {{0.333, 0.333}},
            {{0.50, 1.0e-4}},
            {{1.0e-4, 0.50}},
            {{0.80, 1.0e-4}},
            {{1.0e-4, 0.80}},
            {{0.20, 0.20}},
            {{0.60, 0.20}},
            {{0.20, 0.60}},
            {{0.05, 0.05}},
            {{0.90, 0.01}},
            {{0.95, 1.0e-4}},
            {{0.98, 1.0e-4}},
            {{1.0e-4, 0.95}},
            {{1.0e-4, 0.98}},
            {{0.01, 0.90}},
            {{0.45, 0.45}}}};
        for (const auto &guess : guesses)
        {
            const auto rr = solveRr3_(z, kg, kw, guess[0], guess[1]);
            if (rr.converged)
                return rr;
        }
        return {};
    }

    [[nodiscard]] Rr2 solveRr2_(
        const Composition &z,
        const Composition &K,
        double initialBeta) const
    {
        auto evaluate = [&](double beta, double *dfOut = nullptr) {
            double f = 0.0, df = 0.0;
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double a = K[c] - 1.0;
                const double D = (1.0 - beta) + beta * K[c];
                // Highly asymmetric aqueous/hydrocarbon systems can have
                // equilibrium K values far below 1e-14 for heavy components.
                // At beta=1 the RR denominator then equals K and is still
                // physically positive; rejecting it solely because it is
                // smaller than a generic composition tolerance makes an
                // otherwise well-bracketed G+W flash fail before iteration.
                // Only reject a genuinely non-positive/underflowed denominator.
                if (!(D > 1.0e-300) || !std::isfinite(D))
                    return std::numeric_limits<double>::quiet_NaN();
                f += z[c] * a / D;
                df -= z[c] * a * a / (D * D);
            }
            if (dfOut) *dfOut = df;
            return f;
        };

        const double f0 = evaluate(0.0);
        const double f1 = evaluate(1.0);
        if (!std::isfinite(f0) || !std::isfinite(f1))
            return {};
        if (f0 <= 0.0)
            return {true, 0.0};
        if (f1 >= 0.0)
            return {true, 1.0};

        double lower = 0.0;
        double upper = 1.0;
        double beta = std::clamp(initialBeta, lower, upper);
        for (int iteration = 0; iteration < options_.maximumRachfordRiceIterations; ++iteration)
        {
            double df = 0.0;
            const double f = evaluate(beta, &df);
            if (!std::isfinite(f) || !std::isfinite(df))
                return {};
            if (std::abs(f) < options_.rachfordRiceTolerance)
                return {true, beta};
            if (f > 0.0) lower = beta; else upper = beta;

            double candidate = beta;
            if (std::abs(df) > 1.0e-24)
                candidate = beta - f / df;
            if (!std::isfinite(candidate) || candidate <= lower || candidate >= upper)
                candidate = 0.5 * (lower + upper);
            beta = candidate;
        }
        return {};
    }

    static constexpr std::size_t threePhaseLogKDimension =
        static_cast<std::size_t>(2 * N);

    struct ThreePhaseLogKState
    {
        bool valid{false};
        double norm{std::numeric_limits<double>::infinity()};
        std::array<double, threePhaseLogKDimension> residual{};
        Result result{};
    };

    [[nodiscard]] ThreePhaseLogKState evaluateThreePhaseLogK_(
        double pressure,
        double temperature,
        const Composition &z,
        const std::array<double, threePhaseLogKDimension> &logK,
        double preferredGas,
        double preferredWater) const
    {
        ThreePhaseLogKState state;
        Composition kg{}, kw{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            kg[c] = std::exp(std::clamp(logK[c], -230.0, 230.0));
            kw[c] = std::exp(std::clamp(
                logK[static_cast<std::size_t>(N) + c], -230.0, 230.0));
        }

        const Rr3 rr = solveRr3Robust_(
            z, kg, kw, preferredGas, preferredWater);
        if (!rr.converged)
            return state;
        const double bg = rr.betaGas;
        const double bw = rr.betaWater;
        const double bo = 1.0 - bg - bw;
        if (!(bo >= -options_.phaseFractionTolerance) ||
            !(bg >= -options_.phaseFractionTolerance) ||
            !(bw >= -options_.phaseFractionTolerance))
            return state;

        Composition xo{}, yg{}, xw{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double D = bo + bg * kg[c] + bw * kw[c];
            if (!(D > 0.0) || !std::isfinite(D))
                return state;
            xo[c] = z[c] / D;
            yg[c] = kg[c] * xo[c];
            xw[c] = kw[c] * xo[c];
        }
        normalize_(xo);
        normalize_(yg);
        normalize_(xw);

        const auto to = phaseResult_(
            CompositionalPhase::Oil, pressure, temperature, xo);
        const auto tg = phaseResult_(
            CompositionalPhase::Gas, pressure, temperature, yg);
        const auto tw = phaseResult_(
            CompositionalPhase::Water, pressure, temperature, xw);

        state.norm = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double phiO = std::max(to.fugacityCoefficient[c], 1.0e-300);
            const double phiG = std::max(tg.fugacityCoefficient[c], 1.0e-300);
            const double phiW = std::max(tw.fugacityCoefficient[c], 1.0e-300);
            state.residual[c] = logK[c] - std::log(phiO / phiG);
            state.residual[static_cast<std::size_t>(N) + c] =
                logK[static_cast<std::size_t>(N) + c] - std::log(phiO / phiW);
            state.norm = std::max({
                state.norm,
                std::abs(state.residual[c]),
                std::abs(state.residual[static_cast<std::size_t>(N) + c])});
        }

        state.result.presence = PhasePresence::all();
        state.result.phaseMoleFraction = {bo, bg, bw};
        state.result.composition = {xo, yg, xw};
        state.result.compressibility = {
            to.compressibility, tg.compressibility, tw.compressibility};
        state.result.vaporOilK = kg;
        state.result.waterOilK = kw;
        state.valid = std::isfinite(state.norm);
        return state;
    }

    [[nodiscard]] Result flashThree_(
        double pressure,
        double temperature,
        const Composition &z,
        Composition kg,
        Composition kw,
        double initialBetaGas = 0.15,
        double initialBetaWater = 0.15) const
    {
        Result result;
        result.presence = PhasePresence::all();
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            kg[c] = std::clamp(kg[c], 1.0e-100, 1.0e100);
            kw[c] = std::clamp(kw[c], 1.0e-100, 1.0e100);
        }

        double bg = std::clamp(initialBetaGas, 1.0e-12, 0.999);
        double bw = std::clamp(initialBetaWater, 1.0e-12, 0.999);
        if (bg + bw >= 0.999)
        {
            const double scale = 0.999 / (bg + bw);
            bg *= scale;
            bw *= scale;
        }
        const int ssiSteps = std::min(
            options_.maximumIterations,
            std::max(options_.threePhaseSsiIterationsBeforeNewton, 0));
        int totalIterations = 0;

        // 数值：先用阻尼 SSI 把 K 值推入正确吸引域，再切换到 log-K Newton 加速局部收敛。
        auto oneSsiStep = [&]() -> bool {
            const Rr3 rr = solveRr3Robust_(z, kg, kw, bg, bw);
            if (!rr.converged)
                return false;
            bg = rr.betaGas;
            bw = rr.betaWater;
            const double bo = 1.0 - bg - bw;

            Composition xo{}, yg{}, xw{};
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double D = bo + bg * kg[c] + bw * kw[c];
                if (!(D > 0.0) || !std::isfinite(D))
                    return false;
                xo[c] = z[c] / D;
                yg[c] = kg[c] * xo[c];
                xw[c] = kw[c] * xo[c];
            }
            normalize_(xo);
            normalize_(yg);
            normalize_(xw);

            const auto to = phaseResult_(
                CompositionalPhase::Oil, pressure, temperature, xo);
            const auto tg = phaseResult_(
                CompositionalPhase::Gas, pressure, temperature, yg);
            const auto tw = phaseResult_(
                CompositionalPhase::Water, pressure, temperature, xw);

            double maxResidual = 0.0;
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double targetKg =
                    std::max(to.fugacityCoefficient[c], 1.0e-300) /
                    std::max(tg.fugacityCoefficient[c], 1.0e-300);
                const double targetKw =
                    std::max(to.fugacityCoefficient[c], 1.0e-300) /
                    std::max(tw.fugacityCoefficient[c], 1.0e-300);
                const double dg = std::clamp(
                    std::log(targetKg / kg[c]),
                    -options_.logKStepLimit, options_.logKStepLimit);
                const double dw = std::clamp(
                    std::log(targetKw / kw[c]),
                    -options_.logKStepLimit, options_.logKStepLimit);
                maxResidual = std::max({maxResidual, std::abs(dg), std::abs(dw)});
                kg[c] *= std::exp(0.6 * dg);
                kw[c] *= std::exp(0.6 * dw);
                kg[c] = std::clamp(kg[c], 1.0e-100, 1.0e100);
                kw[c] = std::clamp(kw[c], 1.0e-100, 1.0e100);
            }

            ++totalIterations;
            result.iterations = totalIterations;
            result.phaseMoleFraction = {bo, bg, bw};
            result.composition = {xo, yg, xw};
            result.compressibility = {
                to.compressibility, tg.compressibility, tw.compressibility};
            result.vaporOilK = kg;
            result.waterOilK = kw;

            if (maxResidual < options_.fugacityTolerance)
            {
                result.converged = true;
                result.presence = presenceFromFractions_(result.phaseMoleFraction);
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

        std::array<double, threePhaseLogKDimension> logK{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            logK[c] = std::log(kg[c]);
            logK[static_cast<std::size_t>(N) + c] = std::log(kw[c]);
        }

        for (int iteration = 0;
             iteration < options_.maximumThreePhaseNewtonIterations;
             ++iteration)
        {
            const auto base = evaluateThreePhaseLogK_(
                pressure, temperature, z, logK, bg, bw);
            if (!base.valid)
                break;
            result = base.result;
            result.iterations = ++totalIterations;
            bg = result.phaseMoleFraction[1];
            bw = result.phaseMoleFraction[2];
            if (base.norm < options_.fugacityTolerance)
            {
                result.converged = true;
                result.presence = presenceFromFractions_(result.phaseMoleFraction);
                updateSaturationsFromMoles(pressure, temperature, result);
                return result;
            }

            std::array<
                std::array<double, threePhaseLogKDimension>,
                threePhaseLogKDimension> jacobian{};
            bool jacobianValid = true;
            for (std::size_t column = 0;
                 column < threePhaseLogKDimension;
                 ++column)
            {
                const double step = options_.numericalJacobianStep *
                    std::max(1.0, std::abs(logK[column]));

                auto forwardLogK = logK;
                forwardLogK[column] += step;
                const auto forward = evaluateThreePhaseLogK_(
                    pressure, temperature, z, forwardLogK, bg, bw);
                if (forward.valid)
                {
                    for (std::size_t row = 0;
                         row < threePhaseLogKDimension;
                         ++row)
                    {
                        jacobian[row][column] =
                            (forward.residual[row] - base.residual[row]) / step;
                    }
                    continue;
                }

                // 相边界处正向扰动可能离开可行域。优先改用同尺度后向差分；
                // 两侧都不可评价时放弃本轮 Newton，让已有 SSI globalization 接管。
                auto backwardLogK = logK;
                backwardLogK[column] -= step;
                const auto backward = evaluateThreePhaseLogK_(
                    pressure, temperature, z, backwardLogK, bg, bw);
                if (!backward.valid)
                {
                    jacobianValid = false;
                    break;
                }
                for (std::size_t row = 0;
                     row < threePhaseLogKDimension;
                     ++row)
                {
                    jacobian[row][column] =
                        (base.residual[row] - backward.residual[row]) / step;
                }
            }
            if (!jacobianValid)
                break;

            std::array<double, threePhaseLogKDimension> rhs{};
            std::array<double, threePhaseLogKDimension> delta{};
            for (std::size_t i = 0; i < threePhaseLogKDimension; ++i)
                rhs[i] = -base.residual[i];
            if (!solveDenseLinearSystem_<threePhaseLogKDimension>(
                    jacobian, rhs, delta))
                break;

            bool accepted = false;
            double damping = 1.0;
            for (int lineSearch = 0; lineSearch < 16; ++lineSearch)
            {
                auto trialLogK = logK;
                for (std::size_t i = 0; i < threePhaseLogKDimension; ++i)
                    trialLogK[i] += damping * delta[i];
                const auto trial = evaluateThreePhaseLogK_(
                    pressure, temperature, z, trialLogK, bg, bw);
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

        // 数值：相边界附近数值 Jacobian 可能退化；此时用剩余 SSI 迭代继续全局化，
        // 不因一次 Newton 奇异就立即判定 flash 失败。
        for (int iteration = ssiSteps;
             iteration < options_.maximumIterations;
             ++iteration)
        {
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                kg[c] = std::exp(std::clamp(logK[c], -230.0, 230.0));
                kw[c] = std::exp(std::clamp(
                    logK[static_cast<std::size_t>(N) + c], -230.0, 230.0));
            }
            if (!oneSsiStep())
                return result;
            if (result.converged)
                return result;
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                logK[c] = std::log(kg[c]);
                logK[static_cast<std::size_t>(N) + c] = std::log(kw[c]);
            }
        }
        return result;
    }

    [[nodiscard]] Composition pairInitialK_(
        CompositionalPhase first,
        CompositionalPhase second,
        double pressure,
        double temperature,
        const Composition &z) const
    {
        // SW 的水富集相使用不同的 H2O–溶质作用参数；固定水/油 K 初猜可能远离 G+W/O+W
        // 的正确平衡吸引域，因此优先由当前组成和相角色构造更接近物理状态的初值
        // with many pseudo-components.  Build a phase-aware fugacity-ratio seed
        // whenever SW and a water-rich phase are involved.  Ordinary PR keeps
        // the historical Wilson/fixed-water initialization exactly unchanged.
        if (eos_.usesSoreideWhitson() &&
            (first == CompositionalPhase::Water ||
             second == CompositionalPhase::Water))
        {
            const int water = options_.waterComponent >= 0
                ? options_.waterComponent
                : N - 1;

            auto seedFor = [&](CompositionalPhase phase) {
                Composition seed = z;
                if (phase == CompositionalPhase::Water)
                {
                    // Water-rich trial: retain the feed ordering for dissolved
                    // species but make H2O dominant.  A finite dissolved tail
                    // avoids evaluating log/fugacity derivatives at exact zero.
                    double nonWater = 0.0;
                    for (int i = 0; i < N; ++i)
                    {
                        if (i == water) continue;
                        const std::size_t c = static_cast<std::size_t>(i);
                        seed[c] = std::max(z[c], options_.compositionFloor);
                        nonWater += seed[c];
                    }
                    const double dissolved = 0.02;
                    if (nonWater > 0.0)
                        for (int i = 0; i < N; ++i)
                            if (i != water)
                                seed[static_cast<std::size_t>(i)] *=
                                    dissolved / nonWater;
                    seed[static_cast<std::size_t>(water)] = 1.0 - dissolved;
                }
                else
                {
                    // Hydrocarbon-rich trial: preserve all non-water feed ratios
                    // and use a small but finite H2O content.  This works for
                    // both vapor and liquid hydrocarbon identities.
                    const double waterSeed = 0.01;
                    double nonWater = 0.0;
                    for (int i = 0; i < N; ++i)
                        if (i != water)
                            nonWater += std::max(
                                z[static_cast<std::size_t>(i)],
                                options_.compositionFloor);
                    for (int i = 0; i < N; ++i)
                    {
                        const std::size_t c = static_cast<std::size_t>(i);
                        if (i == water)
                            seed[c] = waterSeed;
                        else
                            seed[c] = nonWater > 0.0
                                ? (1.0 - waterSeed) *
                                      std::max(z[c], options_.compositionFloor) /
                                      nonWater
                                : (1.0 - waterSeed) /
                                      static_cast<double>(std::max(1, N - 1));
                    }
                }
                normalize_(seed);
                return seed;
            };

            const Composition firstSeed = seedFor(first);
            const Composition secondSeed = seedFor(second);
            const auto firstThermo =
                phaseResult_(first, pressure, temperature, firstSeed);
            const auto secondThermo =
                phaseResult_(second, pressure, temperature, secondSeed);

            Composition K{};
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                K[c] = std::clamp(
                    std::max(firstThermo.fugacityCoefficient[c], 1.0e-300) /
                        std::max(secondThermo.fugacityCoefficient[c], 1.0e-300),
                    1.0e-100,
                    1.0e100);
            }
            return K;
        }

        const Composition kg = wilsonGasOilK_(pressure, temperature);
        const Composition kw = initialWaterOilK_(z);
        Composition K{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            auto rel = [&](CompositionalPhase p) {
                if (p == CompositionalPhase::Oil) return 1.0;
                if (p == CompositionalPhase::Gas) return kg[c];
                return kw[c];
            };
            K[c] = rel(second) / rel(first);
        }
        return K;
    }

