/**
 * @file three_phase_flash_allocation.inc
 * @brief 三相 Flash 逐组分守恒分配变量、Gibbs 预处理和三相后备求解实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    static constexpr std::size_t allocationDimension =
        static_cast<std::size_t>(2 * N);
    using AllocationVector = std::array<double, allocationDimension>;
    using AllocationMatrix = std::array<AllocationVector, allocationDimension>;

    struct AllocationState
    {
        bool valid{false};
        double gibbs{0.0};
        double closure{std::numeric_limits<double>::infinity()};
        double residualNorm2{std::numeric_limits<double>::infinity()};
        AllocationVector residual{};
        AllocationVector gradient{};
        std::array<Composition, 3> composition{};
        std::array<double, 3> beta{};
        std::array<double, 3> compressibility{};
    };

    /**
     * @brief 由逐组分软分配变量构造严格守恒的三相候选状态。
     *
     * 数学：对每个组分 i 用 softmax 得到 `q_o+q_g+q_w=1`，并令
     * `n_{p,i}=z_i q_{p,i}`；因此组分总量守恒在参数化层面精确满足，优化只需处理化学势闭合。
     */
    [[nodiscard]] AllocationState evaluateAllocation_(
        double pressure,
        double temperature,
        const Composition &z,
        const AllocationVector &u) const
    {
        AllocationState state;
        std::array<Composition, 3> amount{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double ug = u[c];
            const double uw = u[static_cast<std::size_t>(N) + c];
            const double shift = std::max({0.0, ug, uw});
            const double qo0 = std::exp(-shift);
            const double qg0 = std::exp(ug - shift);
            const double qw0 = std::exp(uw - shift);
            const double qsum = qo0 + qg0 + qw0;
            const double qo = qo0 / qsum;
            const double qg = qg0 / qsum;
            const double qw = qw0 / qsum;
            amount[0][c] = z[c] * qo;
            amount[1][c] = z[c] * qg;
            amount[2][c] = z[c] * qw;
            state.beta[0] += amount[0][c];
            state.beta[1] += amount[1][c];
            state.beta[2] += amount[2][c];
        }

        for (std::size_t p = 0; p < 3; ++p)
        {
            if (!(state.beta[p] > 1.0e-14))
                return state;
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                state.composition[p][c] = amount[p][c] / state.beta[p];
            }
        }

        std::array<Composition, 3> mu{};
        for (std::size_t p = 0; p < 3; ++p)
        {
            const auto thermo = phaseResult_(
                phases_[p], pressure, temperature, state.composition[p]);
            state.compressibility[p] = thermo.compressibility;
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                const double x = std::max(
                    state.composition[p][c], options_.compositionFloor);
                const double phi = std::max(
                    thermo.fugacityCoefficient[c], 1.0e-300);
                mu[p][c] = std::log(x) + std::log(phi);
                state.gibbs += amount[p][c] * mu[p][c];
            }
        }

        state.closure = 0.0;
        state.residualNorm2 = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const std::size_t wrow = static_cast<std::size_t>(N) + c;
            if (!(z[c] > options_.compositionFloor))
            {
                // 数学：总组成中不存在的组分不携带物质，因此不贡献化学势方程或 Gibbs 梯度；
                // 两个分配残差行保持为零，同时避免后续恢复相份额时出现 0/0。
                state.residual[c] = 0.0;
                state.residual[wrow] = 0.0;
                state.gradient[c] = 0.0;
                state.gradient[wrow] = 0.0;
                continue;
            }
            state.residual[c] = mu[1][c] - mu[0][c];
            state.residual[wrow] = mu[2][c] - mu[0][c];
            state.residualNorm2 +=
                state.residual[c] * state.residual[c] +
                state.residual[wrow] * state.residual[wrow];
            state.closure = std::max({
                state.closure,
                std::abs(mu[0][c] - mu[1][c]),
                std::abs(mu[0][c] - mu[2][c]),
                std::abs(mu[1][c] - mu[2][c])});

            const double qo = amount[0][c] / z[c];
            const double qg = amount[1][c] / z[c];
            const double qw = amount[2][c] / z[c];
            const double average =
                qo * mu[0][c] + qg * mu[1][c] + qw * mu[2][c];
            state.gradient[c] = z[c] * qg * (mu[1][c] - average);
            state.gradient[wrow] = z[c] * qw * (mu[2][c] - average);
        }
        state.valid = std::isfinite(state.gibbs) &&
            std::isfinite(state.residualNorm2) &&
            std::isfinite(state.closure);
        return state;
    }

    /**
     * @brief 先以 Gibbs 能为目标对分配变量做 BFGS 型预条件，把候选状态送入正确平衡盆地。
     * @note 该步骤只负责全局化；最终收敛仍由化学势/逸度相等残差判定。
     */
    void gibbsPreconditionAllocation_(
        double pressure,
        double temperature,
        const Composition &z,
        AllocationVector &u,
        AllocationState &state) const
    {
        AllocationMatrix inverseHessian{};
        for (std::size_t i = 0; i < allocationDimension; ++i)
            inverseHessian[i][i] = 1.0;
        state = evaluateAllocation_(pressure, temperature, z, u);
        if (!state.valid)
            return;

        for (int iteration = 0; iteration < 160; ++iteration)
        {
            AllocationVector direction{};
            for (std::size_t i = 0; i < allocationDimension; ++i)
                for (std::size_t j = 0; j < allocationDimension; ++j)
                    direction[i] -= inverseHessian[i][j] * state.gradient[j];

            double directionalDerivative = 0.0;
            double maxStep = 0.0;
            for (std::size_t i = 0; i < allocationDimension; ++i)
            {
                directionalDerivative += state.gradient[i] * direction[i];
                maxStep = std::max(maxStep, std::abs(direction[i]));
            }
            if (!(directionalDerivative < 0.0) || !std::isfinite(directionalDerivative))
            {
                inverseHessian = {};
                directionalDerivative = 0.0;
                maxStep = 0.0;
                for (std::size_t i = 0; i < allocationDimension; ++i)
                {
                    inverseHessian[i][i] = 1.0;
                    direction[i] = -state.gradient[i];
                    directionalDerivative += state.gradient[i] * direction[i];
                    maxStep = std::max(maxStep, std::abs(direction[i]));
                }
            }
            if (maxStep > 4.0)
                for (double &v : direction)
                    v *= 4.0 / maxStep;

            AllocationVector trialU{};
            AllocationState trial;
            double stepLength = 1.0;
            bool accepted = false;
            for (int ls = 0; ls < 24; ++ls)
            {
                for (std::size_t i = 0; i < allocationDimension; ++i)
                    trialU[i] = std::clamp(
                        u[i] + stepLength * direction[i], -40.0, 40.0);
                trial = evaluateAllocation_(
                    pressure, temperature, z, trialU);
                if (trial.valid &&
                    trial.gibbs <= state.gibbs +
                        1.0e-4 * stepLength * directionalDerivative)
                {
                    accepted = true;
                    break;
                }
                stepLength *= 0.5;
            }
            if (!accepted)
                return;

            AllocationVector s{}, y{}, hy{};
            for (std::size_t i = 0; i < allocationDimension; ++i)
            {
                s[i] = trialU[i] - u[i];
                y[i] = trial.gradient[i] - state.gradient[i];
            }
            double ys = 0.0;
            for (std::size_t i = 0; i < allocationDimension; ++i)
                ys += y[i] * s[i];
            if (ys > 1.0e-12)
            {
                for (std::size_t i = 0; i < allocationDimension; ++i)
                    for (std::size_t j = 0; j < allocationDimension; ++j)
                        hy[i] += inverseHessian[i][j] * y[j];
                double yhy = 0.0;
                for (std::size_t i = 0; i < allocationDimension; ++i)
                    yhy += y[i] * hy[i];
                const double coefficient = (ys + yhy) / (ys * ys);
                AllocationMatrix updated = inverseHessian;
                for (std::size_t i = 0; i < allocationDimension; ++i)
                    for (std::size_t j = 0; j < allocationDimension; ++j)
                        updated[i][j] +=
                            coefficient * s[i] * s[j] -
                            (hy[i] * s[j] + s[i] * hy[j]) / ys;
                inverseHessian = updated;
            }
            else
            {
                inverseHessian = {};
                for (std::size_t i = 0; i < allocationDimension; ++i)
                    inverseHessian[i][i] = 1.0;
            }
            u = trialU;
            state = trial;
        }
    }

    /**
     * @brief 用 LM/阻尼 Newton 将三相化学势差收敛到 flash 逸度容差。
     *
     * 数值：目标残差为 `(mu_g-mu_o, mu_w-mu_o)`；LM 阻尼在近奇异 Jacobian 下提供正则化。
     */
    [[nodiscard]] bool refineAllocationChemicalPotentials_(
        double pressure,
        double temperature,
        const Composition &z,
        AllocationVector &u,
        AllocationState &state) const
    {
        state = evaluateAllocation_(pressure, temperature, z, u);
        if (!state.valid)
            return false;
        double lambda = 1.0e-5;
        for (int iteration = 0; iteration < 160; ++iteration)
        {
            if (state.closure < options_.fugacityTolerance)
                return true;

            AllocationMatrix jacobian{};
            for (std::size_t column = 0; column < allocationDimension; ++column)
            {
                const double h = 2.0e-5 * std::max(1.0, std::abs(u[column]));
                AllocationVector plus = u;
                AllocationVector minus = u;
                plus[column] = std::min(40.0, u[column] + h);
                minus[column] = std::max(-40.0, u[column] - h);
                const auto pstate = evaluateAllocation_(
                    pressure, temperature, z, plus);
                const auto mstate = evaluateAllocation_(
                    pressure, temperature, z, minus);
                if (!pstate.valid || !mstate.valid)
                    return false;
                const double denominator = plus[column] - minus[column];
                for (std::size_t row = 0; row < allocationDimension; ++row)
                    jacobian[row][column] =
                        (pstate.residual[row] - mstate.residual[row]) / denominator;
            }

            bool accepted = false;
            for (int lm = 0; lm < 14 && !accepted; ++lm)
            {
                AllocationMatrix normal{};
                AllocationVector rhs{};
                AllocationVector delta{};
                for (std::size_t a = 0; a < allocationDimension; ++a)
                {
                    for (std::size_t b = 0; b < allocationDimension; ++b)
                    {
                        double value = 0.0;
                        for (std::size_t row = 0; row < allocationDimension; ++row)
                            value += jacobian[row][a] * jacobian[row][b];
                        normal[a][b] = value;
                    }
                    double value = 0.0;
                    for (std::size_t row = 0; row < allocationDimension; ++row)
                        value += jacobian[row][a] * state.residual[row];
                    rhs[a] = -value;
                    normal[a][a] += lambda;
                }
                if (!solveDenseLinearSystem_<allocationDimension>(
                        normal, rhs, delta))
                {
                    lambda *= 10.0;
                    continue;
                }
                double maxDelta = 0.0;
                for (double value : delta)
                    maxDelta = std::max(maxDelta, std::abs(value));
                if (maxDelta > 6.0)
                    for (double &value : delta)
                        value *= 6.0 / maxDelta;

                double damping = 1.0;
                for (int ls = 0; ls < 22; ++ls)
                {
                    AllocationVector trialU = u;
                    for (std::size_t i = 0; i < allocationDimension; ++i)
                        trialU[i] = std::clamp(
                            u[i] + damping * delta[i], -40.0, 40.0);
                    const auto trial = evaluateAllocation_(
                        pressure, temperature, z, trialU);
                    if (trial.valid &&
                        trial.residualNorm2 < state.residualNorm2 *
                            (1.0 - 1.0e-5 * damping))
                    {
                        u = trialU;
                        state = trial;
                        lambda = std::max(1.0e-14, lambda * 0.2);
                        accepted = true;
                        break;
                    }
                    damping *= 0.5;
                }
                if (!accepted)
                    lambda *= 10.0;
            }
            if (!accepted)
                return false;
        }
        return state.closure < options_.fugacityTolerance;
    }

    [[nodiscard]] AllocationVector uniformAllocationSeed_(
        double gasBias,
        double waterBias) const
    {
        AllocationVector u{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            u[c] = gasBias;
            u[static_cast<std::size_t>(N) + c] = waterBias;
        }
        return u;
    }

    [[nodiscard]] bool allocationSeedFromPair_(
        double pressure,
        double temperature,
        const Composition &z,
        PhasePresence allowed,
        AllocationVector &u) const
    {
        Result pair = flashRestrictedRecovery_(
            pressure, temperature, z, allowed);
        if (!pair.converged || pair.presence.count() != 2)
            return false;
        const StabilityResult stability = stabilityTest(
            pressure, temperature, z, pair.presence, pair.composition);
        int missing = -1;
        for (CompositionalPhase phase : phases_)
            if (!pair.presence.contains(phase))
            {
                missing = phaseIndex(phase);
                break;
            }
        if (missing < 0 || !stability.valid ||
            !stability.missingPhaseUnstable[static_cast<std::size_t>(missing)])
            return false;

        const Composition &incipient =
            stability.incipientComposition[static_cast<std::size_t>(missing)];
        double epsilon = 1.0e-3;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            if (incipient[c] > 1.0e-30)
                epsilon = std::min(epsilon, 0.15 * z[c] / incipient[c]);
        }
        epsilon = std::max(epsilon, 1.0e-9);
        std::array<Composition, 3> amount{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            if (!(z[c] > options_.compositionFloor))
                continue;
            amount[static_cast<std::size_t>(missing)][c] =
                epsilon * incipient[c];
            const double remaining =
                z[c] - amount[static_cast<std::size_t>(missing)][c];
            if (!(remaining > 0.0))
                return false;
            double pairAmount = 0.0;
            for (CompositionalPhase phase : phases_)
            {
                const int pi = phaseIndex(phase);
                if (pi == missing)
                    continue;
                const std::size_t slot = static_cast<std::size_t>(pi);
                pairAmount += pair.phaseMoleFraction[slot] *
                    pair.composition[slot][c];
            }
            if (!(pairAmount > 0.0))
                return false;
            for (CompositionalPhase phase : phases_)
            {
                const int pi = phaseIndex(phase);
                if (pi == missing)
                    continue;
                const std::size_t slot = static_cast<std::size_t>(pi);
                amount[slot][c] = pair.phaseMoleFraction[slot] *
                    pair.composition[slot][c] * remaining / pairAmount;
            }
        }
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            if (!(z[c] > options_.compositionFloor))
            {
                u[c] = 0.0;
                u[static_cast<std::size_t>(N) + c] = 0.0;
                continue;
            }
            const double qo = std::max(amount[0][c] / z[c], 1.0e-18);
            const double qg = std::max(amount[1][c] / z[c], 1.0e-18);
            const double qw = std::max(amount[2][c] / z[c], 1.0e-18);
            u[c] = std::log(qg / qo);
            u[static_cast<std::size_t>(N) + c] = std::log(qw / qo);
        }
        return true;
    }


    [[nodiscard]] bool allocationSeedFromRawOg_(
        double pressure,
        double temperature,
        const Composition &z,
        AllocationVector &u) const
    {
        Result pair = flashTwo_(
            pressure, temperature, z,
            CompositionalPhase::Oil, CompositionalPhase::Gas,
            wilsonGasOilK_(pressure, temperature));
        if (!pair.converged || presenceFromFractions_(pair.phaseMoleFraction).count() != 2)
            return false;
        pair.presence = presenceFromFractions_(pair.phaseMoleFraction);
        const StabilityResult stability = stabilityTest(
            pressure, temperature, z, pair.presence, pair.composition);
        const int missing = phaseIndex(CompositionalPhase::Water);
        if (!stability.valid ||
            !stability.missingPhaseUnstable[static_cast<std::size_t>(missing)])
            return false;
        const Composition &incipient =
            stability.incipientComposition[static_cast<std::size_t>(missing)];
        double epsilon = 1.0e-3;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            if (incipient[c] > 1.0e-30)
                epsilon = std::min(epsilon, 0.2 * z[c] / incipient[c]);
        }
        epsilon = std::max(epsilon, 1.0e-9);
        std::array<Composition, 3> amount{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            if (!(z[c] > options_.compositionFloor))
            {
                u[c] = 0.0;
                u[static_cast<std::size_t>(N) + c] = 0.0;
                continue;
            }
            amount[2][c] = epsilon * incipient[c];
            const double remaining = z[c] - amount[2][c];
            const double no = pair.phaseMoleFraction[0] * pair.composition[0][c];
            const double ng = pair.phaseMoleFraction[1] * pair.composition[1][c];
            const double sum = no + ng;
            if (!(remaining > 0.0) || !(sum > 0.0))
                return false;
            amount[0][c] = no * remaining / sum;
            amount[1][c] = ng * remaining / sum;
            const double qo = std::max(amount[0][c] / z[c], 1.0e-18);
            const double qg = std::max(amount[1][c] / z[c], 1.0e-18);
            const double qw = std::max(amount[2][c] / z[c], 1.0e-18);
            u[c] = std::log(qg / qo);
            u[static_cast<std::size_t>(N) + c] = std::log(qw / qo);
        }
        return true;
    }

    [[nodiscard]] bool allocationSeedFromSingle_(
        double pressure,
        double temperature,
        const Composition &z,
        AllocationVector &u) const
    {
        const PhasePresence active = chooseSinglePhase_(pressure, temperature, z);
        const Result single = singlePhaseResult_(pressure, temperature, z, active);
        const StabilityResult stability = stabilityTest(
            pressure, temperature, z, active, single.composition);
        if (!stability.valid || stability.stable)
            return false;
        std::array<double, 3> epsilon{};
        int missingCount = 0;
        for (CompositionalPhase phase : phases_)
        {
            const int pi = phaseIndex(phase);
            if (!active.contains(phase) &&
                stability.missingPhaseUnstable[static_cast<std::size_t>(pi)])
            {
                epsilon[static_cast<std::size_t>(pi)] = 1.0e-4;
                ++missingCount;
            }
        }
        if (missingCount == 0)
            return false;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            for (CompositionalPhase phase : phases_)
            {
                const int pi = phaseIndex(phase);
                const std::size_t slot = static_cast<std::size_t>(pi);
                if (!(epsilon[slot] > 0.0))
                    continue;
                epsilon[slot] = std::min(
                    epsilon[slot],
                    0.08 * z[c] /
                        std::max(stability.incipientComposition[slot][c], 1.0e-30));
            }
        }
        std::array<Composition, 3> amount{};
        const std::size_t activeSlot = static_cast<std::size_t>(
            phaseIndex(firstActive_(active)));
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            if (!(z[c] > options_.compositionFloor))
            {
                u[c] = 0.0;
                u[static_cast<std::size_t>(N) + c] = 0.0;
                continue;
            }
            double removed = 0.0;
            for (CompositionalPhase phase : phases_)
            {
                const std::size_t slot = static_cast<std::size_t>(phaseIndex(phase));
                if (!(epsilon[slot] > 0.0))
                    continue;
                amount[slot][c] = epsilon[slot] *
                    stability.incipientComposition[slot][c];
                removed += amount[slot][c];
            }
            if (!(z[c] > removed))
                return false;
            amount[activeSlot][c] = z[c] - removed;
            const double qo = std::max(amount[0][c] / z[c], 1.0e-18);
            const double qg = std::max(amount[1][c] / z[c], 1.0e-18);
            const double qw = std::max(amount[2][c] / z[c], 1.0e-18);
            u[c] = std::log(qg / qo);
            u[static_cast<std::size_t>(N) + c] = std::log(qw / qo);
        }
        return true;
    }


    [[nodiscard]] AllocationVector thermodynamicAllocationSeed_(
        double pressure,
        double temperature,
        const Composition &z,
        double scale,
        double gasShift,
        double waterShift) const
    {
        const Composition kg = wilsonGasOilK_(pressure, temperature);
        const Composition kw = initialWaterOilK_(z);
        AllocationVector u{};
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            u[c] = std::clamp(
                scale * std::log(std::clamp(kg[c], 1.0e-30, 1.0e30)) +
                    gasShift,
                -30.0, 30.0);
            u[static_cast<std::size_t>(N) + c] = std::clamp(
                scale * std::log(std::clamp(kw[c], 1.0e-30, 1.0e30)) +
                    waterShift,
                -30.0, 30.0);
        }
        return u;
    }

    /**
     * @brief SW 三相分配恢复：从两相/单相/热力学启发式多组种子搜索最低 Gibbs 的三相平衡。
     * @note 仅在常规 flash 失败后使用，不改变 EOS、BIP 或相稳定性判据。
     */
    [[nodiscard]] Result flashThreeAllocationFallback_(
        double pressure,
        double temperature,
        const Composition &z) const
    {
        constexpr std::size_t maximumSeeds = 24;
        std::array<AllocationVector, maximumSeeds> seeds{};
        std::array<bool, maximumSeeds> valid{};
        std::size_t count = 0;
        valid[count] = allocationSeedFromRawOg_(
            pressure, temperature, z, seeds[count]);
        ++count;
        valid[count] = allocationSeedFromSingle_(
            pressure, temperature, z, seeds[count]);
        ++count;
        const std::array<PhasePresence, 3> pairs{{
            PhasePresence(PhasePresence::bit(CompositionalPhase::Oil) |
                          PhasePresence::bit(CompositionalPhase::Gas)),
            PhasePresence(PhasePresence::bit(CompositionalPhase::Oil) |
                          PhasePresence::bit(CompositionalPhase::Water)),
            PhasePresence(PhasePresence::bit(CompositionalPhase::Gas) |
                          PhasePresence::bit(CompositionalPhase::Water))}};
        for (const PhasePresence pair : pairs)
        {
            valid[count] = allocationSeedFromPair_(
                pressure, temperature, z, pair, seeds[count]);
            ++count;
        }
        const std::array<std::array<double, 3>, 6> thermoSeeds{{
            {{1.0, 0.0, 0.0}},
            {{0.5, 0.0, 0.0}},
            {{1.5, 0.0, 0.0}},
            {{1.0, -1.5, -1.5}},
            {{1.0, 1.5, -1.5}},
            {{1.0, -1.5, 1.5}}}};
        for (const auto &spec : thermoSeeds)
        {
            seeds[count] = thermodynamicAllocationSeed_(
                pressure, temperature, z, spec[0], spec[1], spec[2]);
            valid[count] = true;
            ++count;
        }
        const std::array<std::array<double, 2>, 12> generic{{
            {{0.0, 0.0}}, {{1.0, -2.0}}, {{2.0, -4.0}},
            {{2.0, -2.0}}, {{1.0, -4.0}}, {{-1.0, -1.0}},
            {{-2.0, -2.0}}, {{-2.0, 1.0}}, {{-1.0, 2.0}},
            {{2.0, 1.0}}, {{1.0, 2.0}}, {{-2.0, 2.0}}}};
        for (const auto &bias : generic)
        {
            seeds[count] = uniformAllocationSeed_(bias[0], bias[1]);
            valid[count] = true;
            ++count;
        }

        Result best;
        double bestG = std::numeric_limits<double>::infinity();
        for (std::size_t seedIndex = 0; seedIndex < count; ++seedIndex)
        {
            if (!valid[seedIndex])
                continue;
            AllocationVector u = seeds[seedIndex];
            AllocationState state;
            gibbsPreconditionAllocation_(
                pressure, temperature, z, u, state);
            if (!state.valid)
                continue;
            if (!refineAllocationChemicalPotentials_(
                    pressure, temperature, z, u, state))
                continue;
            if (!(state.beta[0] > options_.phaseFractionTolerance) ||
                !(state.beta[1] > options_.phaseFractionTolerance) ||
                !(state.beta[2] > options_.phaseFractionTolerance))
                continue;
            if (!(state.gibbs < bestG))
                continue;

            Result candidate;
            candidate.converged = true;
            candidate.presence = PhasePresence::all();
            candidate.phaseMoleFraction = state.beta;
            candidate.composition = state.composition;
            candidate.compressibility = state.compressibility;
            deriveReferenceK_(candidate);
            updateSaturationsFromMoles(
                pressure, temperature, candidate);
            if (!candidate.converged)
                continue;
            best = candidate;
            bestG = state.gibbs;
        }
        return best;
    }

    /**
     * @brief 在固定活动相集上求平衡，并用缺失相稳定性测试驱动相出现/消失。
     *
     * 状态：相数转换只发生在该层；底层两相/三相求解器只负责给定相集的平衡。
     */
    [[nodiscard]] Result equilibrateActiveSet_(
        double pressure,
        double temperature,
        const Composition &z,
        PhasePresence active,
        Result seed,
        bool reuseConvergedSeed = true) const
    {
        if (active.empty())
            active = chooseSinglePhase_(pressure, temperature, z);

        canonicalize_(pressure, temperature, seed);
        if (seed.presence.bits() != active.bits())
            active = seed.presence;
        Result lastPhysical = seed;

        // 数值：PR/CPA 的收敛液相在标准化后可能换到另一个公开槽位；物理解本身仍有效。
        // 第一次外循环直接复用该收敛状态做缺失相稳定性测试，只有活动相集真正改变后才重新 flash，
        // 避免仅因广义 RR 参考方向变化而制造一次无意义的失败。
        bool useConvergedSeed = reuseConvergedSeed && seed.converged &&
            seed.presence.bits() == active.bits() &&
            !eos_.usesSoreideWhitson();

        for (int outer = 0; outer < 8; ++outer)
        {
            Result current;
            if (useConvergedSeed)
            {
                current = seed;
                useConvergedSeed = false;
            }
            else if (active.count() == 1)
            {
                current = singlePhaseResult_(pressure, temperature, z, active);
            }
            else if (active.count() == 2)
            {
                const auto pair = activePair_(active);
                Composition K = pairInitialK_(
                    pair.first, pair.second, pressure, temperature, z);
                if (seed.converged && seed.presence.bits() == active.bits())
                {
                    const auto &firstComposition = seed.composition[
                        static_cast<std::size_t>(phaseIndex(pair.first))];
                    const auto &secondComposition = seed.composition[
                        static_cast<std::size_t>(phaseIndex(pair.second))];
                    for (int i = 0; i < N; ++i)
                    {
                        const std::size_t c = static_cast<std::size_t>(i);
                        K[c] = std::clamp(
                            secondComposition[c] /
                                std::max(firstComposition[c], options_.compositionFloor),
                            1.0e-100, 1.0e100);
                    }
                }
                current = flashTwo_(
                    pressure, temperature, z, pair.first, pair.second, K);
                if (!current.converged)
                {
                    current = singlePhaseResult_(pressure, temperature, z, chooseSinglePhase_(pressure, temperature, z));
                }
            }
            else
            {
                Composition kg = seed.vaporOilK;
                Composition kw = seed.waterOilK;
                bool kgValid = true, kwValid = true;
                for (int i = 0; i < N; ++i)
                {
                    kgValid = kgValid && kg[static_cast<std::size_t>(i)] > 0.0;
                    kwValid = kwValid && kw[static_cast<std::size_t>(i)] > 0.0;
                }
                if (!kgValid) kg = wilsonGasOilK_(pressure, temperature);
                if (!kwValid) kw = initialWaterOilK_(z);
                current = flashThree_(pressure, temperature, z, kg, kw);
                if (!current.converged)
                {
                    lastPhysical.converged = false;
                    canonicalize_(pressure, temperature, lastPhysical);
                    return lastPhysical;
                }
            }

            PhasePresence fractionPresence = presenceFromFractions_(current.phaseMoleFraction);
            if (fractionPresence.empty())
                fractionPresence = chooseSinglePhase_(pressure, temperature, z);
            if (fractionPresence.bits() != active.bits())
            {
                active = fractionPresence;
                seed = current;
                canonicalize_(pressure, temperature, seed);
                lastPhysical = seed;
                useConvergedSeed = seed.converged &&
                    seed.presence.bits() == active.bits() &&
                    !eos_.usesSoreideWhitson();
                continue;
            }
            canonicalize_(pressure, temperature, current);
            lastPhysical = current;

            if (active.count() == 3)
            {
                canonicalize_(pressure, temperature, current);
                return current;
            }

            const StabilityResult stability = stabilityTest(
                pressure, temperature, z, active, current.composition);
            if (!stability.valid)
            {
                // 数值失败不是热力学“稳定”证书。将本 active-set 路径标为失败，
                // 由上层直接三相/稳定约化子集 fallback 继续搜索，禁止把未知状态
                // 以 converged=true 静默返回。
                current.converged = false;
                canonicalize_(pressure, temperature, current);
                return current;
            }
            if (stability.stable)
            {
                canonicalize_(pressure, temperature, current);
                return current;
            }

            bool added = false;

            // CPA can enter a one-density-root region in which a Water-only
            // equilibrium reports numerically identical missing Oil and Gas
            // TPD directions. Releasing both labels at once creates an
            // artificial W -> O+G+W jump and can miss the physical O+W
            // bifurcation. When both trial strength and incipient composition
            // coincide, treat them as one nonaqueous instability direction:
            // release the liquid-root Oil role first, re-equilibrate, then
            // rerun stability. A genuinely distinct Gas instability will
            // therefore be released on the following outer iteration.
            bool deduplicatedCpaNonaqueousDirection = false;
            if (eos_.usesCubicPlusAssociation() &&
                active.count() == 1 &&
                active.contains(CompositionalPhase::Water) &&
                stability.missingPhaseUnstable[static_cast<std::size_t>(
                    phaseIndex(CompositionalPhase::Oil))] &&
                stability.missingPhaseUnstable[static_cast<std::size_t>(
                    phaseIndex(CompositionalPhase::Gas))])
            {
                const std::size_t oilSlot = static_cast<std::size_t>(
                    phaseIndex(CompositionalPhase::Oil));
                const std::size_t gasSlot = static_cast<std::size_t>(
                    phaseIndex(CompositionalPhase::Gas));
                const double oilTrial = stability.trialSum[oilSlot];
                const double gasTrial = stability.trialSum[gasSlot];
                const double trialScale = std::max({
                    1.0, std::abs(oilTrial), std::abs(gasTrial)});
                const double trialTolerance = std::max(
                    10.0 * options_.stabilityTolerance * trialScale,
                    1.0e-12);

                double incipientL1 = 0.0;
                for (int i = 0; i < N; ++i)
                {
                    const std::size_t component =
                        static_cast<std::size_t>(i);
                    incipientL1 += std::abs(
                        stability.incipientComposition[oilSlot][component] -
                        stability.incipientComposition[gasSlot][component]);
                }
                const double compositionTolerance = std::max(
                    10.0 * options_.coincidentPhaseCompositionTolerance,
                    1.0e-10);

                if (std::isfinite(oilTrial) &&
                    std::isfinite(gasTrial) &&
                    std::abs(oilTrial - gasTrial) <= trialTolerance &&
                    incipientL1 <= compositionTolerance)
                {
                    active.add(CompositionalPhase::Oil);
                    seed.composition[oilSlot] =
                        stability.incipientComposition[oilSlot];
                    added = true;
                    deduplicatedCpaNonaqueousDirection = true;
                }
            }

            if (!deduplicatedCpaNonaqueousDirection)
            {
                for (CompositionalPhase p : phases_)
                {
                    const int pi = phaseIndex(p);
                    if (!active.contains(p) &&
                        stability.missingPhaseUnstable.at(
                            static_cast<std::size_t>(pi)))
                    {
                        active.add(p);
                        seed.composition[static_cast<std::size_t>(pi)] =
                            stability.incipientComposition.at(
                                static_cast<std::size_t>(pi));
                        added = true;
                    }
                }
            }
            if (!added)
            {
                // valid && !stable 必须至少指出一个缺失相不稳定；若没有，说明
                // stability 结果内部不自洽，不能把当前 reduced state 当作物理解。
                current.converged = false;
                canonicalize_(pressure, temperature, current);
                return current;
            }

            deriveReferenceK_(seed);
            seed.presence = active;
            useConvergedSeed = false;
        }
        canonicalize_(pressure, temperature, seed);
        return seed;
    }

