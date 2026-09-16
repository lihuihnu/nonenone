/**
 * @file three_phase_flash_stability.inc
 * @brief 三相 Flash 切平面稳定性多起点、Newton/SSI 和线性求解辅助实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    /**
     * @brief 生成文献常用的多起点组成用于切平面稳定性搜索。
     *
     * Pang (2017), following Li & Firoozabadi, recommends Nc+4 stability
     * initializations: Wilson K, reciprocal Wilson K, cube-root Wilson K,
     * reciprocal cube-root Wilson K, and Nc component-enriched trials.  The
     * trial equations below only need compositions, so we form Y_i=K_i z_i
     * and normalize.  Two inexpensive supplementary seeds (z and the active
     * reference phase) are retained because they are useful close to a
     * critical/trivial stationary point.
     */
    struct StabilitySeeds_ final
    {
        std::array<Composition, static_cast<std::size_t>(N + 7)> values{};
        std::size_t size{0};

        void push(const Composition &composition)
        {
            if (size >= values.size())
                throw std::logic_error("Three-phase stability seed capacity exceeded.");
            values[size++] = composition;
        }
    };

    [[nodiscard]] StabilitySeeds_ stabilitySeeds_(
        CompositionalPhase candidate,
        double pressure,
        double temperature,
        const Composition &z,
        const Composition &referenceComposition) const
    {
        (void)candidate;
        StabilitySeeds_ seeds;

        const Composition wilson = wilsonGasOilK_(pressure, temperature);
        auto appendKSeed = [&](double exponent) {
            Composition seed{};
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double k = std::clamp(wilson[c], 1.0e-100, 1.0e100);
                seed[c] = z[c] * std::pow(k, exponent);
            }
            normalize_(seed);
            seeds.push(seed);
        };

        appendKSeed(1.0);
        appendKSeed(-1.0);
        appendKSeed(1.0 / 3.0);
        appendKSeed(-1.0 / 3.0);

        if constexpr (N == 1)
        {
            seeds.push(Composition{1.0});
        }
        else
        {
            for (int enriched = 0; enriched < N; ++enriched)
            {
                Composition seed{};
                const double remainder = 0.1 / static_cast<double>(N - 1);
                seed.fill(remainder);
                seed[static_cast<std::size_t>(enriched)] = 0.9;
                seeds.push(seed);
            }
        }

        Composition zSeed = z;
        Composition referenceSeed = referenceComposition;
        normalize_(zSeed);
        normalize_(referenceSeed);
        seeds.push(zSeed);
        seeds.push(referenceSeed);

        // Water-containing systems benefit from an almost-pure H2O trial; this
        // is the water-specific extension proposed in Pang's later chapter and
        // costs only one additional TPD search.
        if (options_.waterComponent >= 0 && N > 1)
        {
            Composition waterSeed{};
            const std::size_t w = static_cast<std::size_t>(options_.waterComponent);
            const double trace = 1.0e-6 / static_cast<double>(N - 1);
            waterSeed.fill(trace);
            waterSeed[w] = 1.0 - 1.0e-6;
            seeds.push(waterSeed);
        }
        return seeds;
    }

    template <std::size_t M>
    [[nodiscard]] static bool solveDenseLinearSystem_(
        std::array<std::array<double, M>, M> matrix,
        std::array<double, M> rhs,
        std::array<double, M> &solution)
    {
        solution.fill(0.0);
        for (std::size_t column = 0; column < M; ++column)
        {
            std::size_t pivot = column;
            for (std::size_t row = column + 1; row < M; ++row)
            {
                if (std::abs(matrix[row][column]) >
                    std::abs(matrix[pivot][column]))
                    pivot = row;
            }
            if (!(std::abs(matrix[pivot][column]) > 1.0e-13))
                return false;
            if (pivot != column)
            {
                std::swap(matrix[pivot], matrix[column]);
                std::swap(rhs[pivot], rhs[column]);
            }

            for (std::size_t row = column + 1; row < M; ++row)
            {
                const double factor = matrix[row][column] / matrix[column][column];
                if (!std::isfinite(factor))
                    return false;
                for (std::size_t j = column; j < M; ++j)
                    matrix[row][j] -= factor * matrix[column][j];
                rhs[row] -= factor * rhs[column];
            }
        }

        for (std::size_t reverse = 0; reverse < M; ++reverse)
        {
            const std::size_t row = M - 1 - reverse;
            double value = rhs[row];
            for (std::size_t j = row + 1; j < M; ++j)
                value -= matrix[row][j] * solution[j];
            if (!(std::abs(matrix[row][row]) > 1.0e-13))
                return false;
            solution[row] = value / matrix[row][row];
            if (!std::isfinite(solution[row]))
                return false;
        }
        return true;
    }

    struct StabilityResidual
    {
        bool valid{false};
        double norm{std::numeric_limits<double>::infinity()};
        double sum{0.0};
        Composition composition{};
        Composition residual{};
    };

    [[nodiscard]] StabilityResidual evaluateStabilityResidual_(
        CompositionalPhase candidate,
        double pressure,
        double temperature,
        const std::array<double, N> &logReference,
        const Composition &logW) const
    {
        StabilityResidual state;
        const double maxLogW = *std::max_element(logW.begin(), logW.end());
        if (!std::isfinite(maxLogW))
            return state;

        double scaledSum = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double value = std::exp(logW[c] - maxLogW);
            if (!(value > 0.0) || !std::isfinite(value))
                return state;
            state.composition[c] = value;
            scaledSum += value;
        }
        if (!(scaledSum > 0.0) || !std::isfinite(scaledSum))
            return state;
        for (double &value : state.composition)
            value /= scaledSum;

        const double logSum = maxLogW + std::log(scaledSum);
        if (logSum > std::log(std::numeric_limits<double>::max()))
            state.sum = std::numeric_limits<double>::infinity();
        else
            state.sum = std::exp(logSum);

        const auto thermo = phaseResult_(
            candidate, pressure, temperature, state.composition);
        state.norm = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double phi = thermo.fugacityCoefficient[c];
            if (!(phi > 0.0) || !std::isfinite(phi))
                return StabilityResidual{};
            state.residual[c] = logW[c] + std::log(phi) - logReference[c];
            state.norm = std::max(state.norm, std::abs(state.residual[c]));
        }
        state.valid = std::isfinite(state.norm);
        return state;
    }

    [[nodiscard]] Trial finalizeStabilityTrial_(
        const StabilityResidual &state,
        const Composition &referenceComposition) const
    {
        if (!state.valid)
            return {};
        double referenceDistance = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            referenceDistance +=
                std::abs(state.composition[c] - referenceComposition[c]);
        }
        const bool trivial =
            referenceDistance < 1.0e-7 &&
            std::abs(state.sum - 1.0) < 1.0e-7;
        return {true, true, trivial, state.sum, state.composition};
    }

    [[nodiscard]] Trial stabilityTrial_(
        CompositionalPhase candidate,
        double pressure,
        double temperature,
        const std::array<double, N> &logReference,
        const Composition &referenceComposition,
        Composition seed) const
    {
        normalize_(seed);
        Composition logW{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            logW[c] = std::log(std::max(seed[c], options_.compositionFloor));
        }

        // A few inexpensive SSI steps place the iterate in the attraction
        // basin. High-pressure water/oil/gas TPD problems are well known to
        // converge very slowly with SSI alone, so a numerical Newton solve is
        // used as the robust second-order stage.
        const int ssiSteps = std::min(
            options_.maximumStabilityIterations,
            std::max(options_.stabilitySsiIterationsBeforeNewton, 0));
        for (int iteration = 0; iteration < ssiSteps; ++iteration)
        {
            const auto current = evaluateStabilityResidual_(
                candidate, pressure, temperature, logReference, logW);
            if (!current.valid)
                return {};
            if (current.norm < options_.stabilityTolerance)
                return finalizeStabilityTrial_(current, referenceComposition);

            const auto thermo = phaseResult_(
                candidate, pressure, temperature, current.composition);
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double phi = std::max(
                    thermo.fugacityCoefficient[c], 1.0e-300);
                logW[c] = logReference[c] - std::log(phi);
            }
        }

        for (int iteration = 0;
             iteration < options_.maximumStabilityNewtonIterations;
             ++iteration)
        {
            const auto base = evaluateStabilityResidual_(
                candidate, pressure, temperature, logReference, logW);
            if (!base.valid)
                return {};
            if (base.norm < options_.stabilityTolerance)
                return finalizeStabilityTrial_(base, referenceComposition);

            std::array<std::array<double, N>, N> jacobian{};
            for (int column = 0; column < N; ++column)
            {
                const std::size_t j = static_cast<std::size_t>(column);
                Composition perturbed = logW;
                const double step = options_.numericalJacobianStep *
                    std::max(1.0, std::abs(logW[j]));
                perturbed[j] += step;
                const auto shifted = evaluateStabilityResidual_(
                    candidate, pressure, temperature, logReference, perturbed);
                if (!shifted.valid)
                    return {};
                for (int row = 0; row < N; ++row)
                {
                    const std::size_t i = static_cast<std::size_t>(row);
                    jacobian[i][j] =
                        (shifted.residual[i] - base.residual[i]) / step;
                }
            }

            Composition rhs{};
            Composition delta{};
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
                Composition candidateLogW = logW;
                for (int i = 0; i < N; ++i)
                {
                    const std::size_t c = static_cast<std::size_t>(i);
                    candidateLogW[c] += damping * delta[c];
                }
                const auto trial = evaluateStabilityResidual_(
                    candidate, pressure, temperature, logReference,
                    candidateLogW);
                if (trial.valid && trial.norm < base.norm)
                {
                    logW = candidateLogW;
                    accepted = true;
                    break;
                }
                damping *= 0.5;
            }
            if (!accepted)
                break;
        }

        // Fall back to the remaining SSI budget rather than returning a false
        // stable result when Newton is singular near a critical point.
        const int remainingSsi = std::max(
            0, options_.maximumStabilityIterations - ssiSteps);
        for (int iteration = 0; iteration < remainingSsi; ++iteration)
        {
            const auto current = evaluateStabilityResidual_(
                candidate, pressure, temperature, logReference, logW);
            if (!current.valid)
                return {};
            if (current.norm < options_.stabilityTolerance)
                return finalizeStabilityTrial_(current, referenceComposition);
            const auto thermo = phaseResult_(
                candidate, pressure, temperature, current.composition);
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                logW[c] = logReference[c] - std::log(std::max(
                    thermo.fugacityCoefficient[c], 1.0e-300));
            }
        }
        return {};
    }

