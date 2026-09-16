/**
 * @file three_phase_flash_canonicalization.inc
 * @brief 三相 Flash O/G/W 相角色规范化、物理证书和最终一致性实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    void canonicalize_(double pressure, double temperature, Result &result) const
    {
        PhasePresence p = presenceFromFractions_(result.phaseMoleFraction);
        if (p.empty())
            p = result.presence.empty() ? PhasePresence::oilOnly() : result.presence;
        result.presence = p;
        double sum = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            const auto phase = phases_[static_cast<std::size_t>(i)];
            double &beta = result.phaseMoleFraction[static_cast<std::size_t>(i)];
            if (!p.contains(phase) || beta <= options_.phaseFractionTolerance || !std::isfinite(beta))
            {
                beta = 0.0;
                result.presence.remove(phase);
            }
            else
            {
                sum += beta;
            }
        }
        if (!(sum > 0.0))
        {
            result.presence = chooseSinglePhase_(pressure, temperature, result.composition[0]);
            result.phaseMoleFraction = {0.0, 0.0, 0.0};
            result.phaseMoleFraction[static_cast<std::size_t>(phaseIndex(firstActive_(result.presence)))] = 1.0;
        }
        else
        {
            for (double &beta : result.phaseMoleFraction)
                beta /= sum;
        }

        if (!eos_.usesSoreideWhitson())
            canonicalizeLiquidRoles_(result);
        result.presence = presenceFromFractions_(result.phaseMoleFraction);
        refreshCanonicalThermodynamics_(pressure, temperature, result);
    }

    void refreshCanonicalThermodynamics_(
        double pressure,
        double temperature,
        Result &result) const
    {
        for (CompositionalPhase phase : phases_)
        {
            const int pi = phaseIndex(phase);
            if (!result.presence.contains(phase))
                continue;
            normalize_(result.composition[static_cast<std::size_t>(pi)]);
            const auto thermo = phaseResult_(
                phase, pressure, temperature,
                result.composition[static_cast<std::size_t>(pi)]);
            result.compressibility[static_cast<std::size_t>(pi)] = thermo.compressibility;
        }
        deriveReferenceK_(result);
        updateSaturationsFromMoles(pressure, temperature, result);
    }

    [[nodiscard]] bool waterIsDominant_(const Composition &composition) const noexcept
    {
        if (options_.waterComponent < 0)
            return false;
        // A public Water role also selects the configured aqueous density and
        // viscosity closures.  Composition enrichment alone is therefore not
        // sufficient: an H2O-largest hydrocarbon-rich PR phase must remain a
        // cubic O/G phase when it lies outside that explicit aqueous domain.
        if (!eos_.aqueousVolumeCompositionSupported(composition))
            return false;
        const std::size_t w = static_cast<std::size_t>(options_.waterComponent);
        double largestOther = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            if (c != w)
                largestOther = std::max(largestOther, composition[c]);
        }
        return composition[w] >= largestOther;
    }

    /**
     * @brief 按 H2O 富集程度规范化两个 PR 液相的物理身份。
     *
     * Oil and water-rich phases both use the PR liquid root; their identity is
     * therefore compositional, not a different EOS root.  A restricted G+W
     * solve can numerically converge to a hydrocarbon-rich liquid in the W
     * slot.  Without relabeling, a subsequent "missing oil" stability trial
     * sees that same liquid as a trivial solution and can suppress the true
     * water-rich third phase.  When both liquids exist, the one with more H2O
     * is W.  With only one liquid, H2O must be its dominant component to label
     * it W; otherwise it is O.
     */
    void canonicalizeLiquidRoles_(Result &result) const
    {
        if (options_.waterComponent < 0)
            return;

        constexpr std::size_t oil = 0;
        constexpr std::size_t water = 2;
        const bool hasOil =
            result.phaseMoleFraction[oil] > options_.phaseFractionTolerance;
        const bool hasWater =
            result.phaseMoleFraction[water] > options_.phaseFractionTolerance;

        bool swapLiquidSlots = false;
        if (hasOil && hasWater)
        {
            const std::size_t w = static_cast<std::size_t>(options_.waterComponent);
            swapLiquidSlots =
                result.composition[oil][w] > result.composition[water][w];
        }
        else if (hasOil != hasWater)
        {
            const std::size_t active = hasOil ? oil : water;
            const bool shouldBeWater = waterIsDominant_(result.composition[active]);
            swapLiquidSlots = (hasOil && shouldBeWater) || (hasWater && !shouldBeWater);
        }

        if (!swapLiquidSlots)
            return;

        std::swap(result.phaseMoleFraction[oil], result.phaseMoleFraction[water]);
        std::swap(result.saturation[oil], result.saturation[water]);
        std::swap(result.composition[oil], result.composition[water]);
        std::swap(result.compressibility[oil], result.compressibility[water]);
        std::swap(result.molarDensity[oil], result.molarDensity[water]);
    }

    /**
     * @brief 按 Wilson 挥发性规范化共存非水相身份。
     *
     * Oil and gas slots can exchange numerical identities when both phase
     * compositions lie in a one-root or near-critical region.  The converged
     * equilibrium itself is still valid, but fixed slot labels then invert
     * experimental oil/gas comparisons and can confuse downstream phase-state
     * logic.  Use the same Wilson volatility ordering already employed by the
     * flash stability seeds: the gas composition must be enriched, in aggregate,
     * in components with larger Wilson K values.  This remains meaningful when
     * liquid/vapour cubic roots merge, unlike a root-index label.  Near zero
     * volatility contrast the phase identity is ill-conditioned, so retain the
     * existing ordering instead of introducing label chatter.
     */
    void swapPhaseSlots_(Result &result, std::size_t a, std::size_t b) const
    {
        if (a == b)
            return;
        std::swap(result.phaseMoleFraction[a], result.phaseMoleFraction[b]);
        std::swap(result.saturation[a], result.saturation[b]);
        std::swap(result.composition[a], result.composition[b]);
        std::swap(result.compressibility[a], result.compressibility[b]);
        std::swap(result.molarDensity[a], result.molarDensity[b]);
    }

    /**
     * @brief 将物理上最富水的相映射到专用 W 输出槽位。
     *
     * The old liquid-only rule compared O and W because both are liquid roots.
     * A fully compositional three-phase solve can, however, return the water-rich
     * physical phase in the numerical G slot as well.  That is harmless only for
     * a completely phase-symmetric EOS; it changes the cubic root for PR/CPA and
     * changes the aqueous BIP model for Soreide-Whitson.  At the public flash
     * boundary therefore compare every active slot by H2O enrichment.  If W is
     * present, the most H2O-rich active composition owns the W identity.  If W
     * is absent, create the W identity only when the richest active composition
     * is itself water-dominant, preserving ordinary O/G states with trace water.
     *
     * @return Source slot moved into W, or -1 when no water-role change occurred.
     */
    [[nodiscard]] int canonicalizeGlobalWaterRole_(Result &result) const
    {
        if (options_.waterComponent < 0)
            return -1;

        constexpr std::size_t waterSlot = 2;
        const std::size_t wc = static_cast<std::size_t>(options_.waterComponent);
        std::size_t richest = waterSlot;
        double richestWater = -1.0;
        int activeCount = 0;
        bool hasWaterSlot = false;

        for (std::size_t slot = 0; slot < 3; ++slot)
        {
            if (!(result.phaseMoleFraction[slot] > options_.phaseFractionTolerance))
                continue;
            ++activeCount;
            if (slot == waterSlot)
                hasWaterSlot = true;
            const double xw = result.composition[slot][wc];
            if (xw > richestWater)
            {
                richestWater = xw;
                richest = slot;
            }
        }

        if (activeCount == 0 || richest == waterSlot)
            return -1;

        if (!hasWaterSlot && !waterIsDominant_(result.composition[richest]))
            return -1;

        const double currentWater = hasWaterSlot
            ? result.composition[waterSlot][wc]
            : -1.0;
        constexpr double relativeWaterRoleTolerance = 1.0e-10;
        if (hasWaterSlot &&
            richestWater <= currentWater +
                relativeWaterRoleTolerance * std::max(1.0, std::abs(richestWater)))
        {
            return -1;
        }

        swapPhaseSlots_(result, richest, waterSlot);
        return static_cast<int>(richest);
    }

    /** @brief 统一验证任何对外 `converged=true` flash 的物料与逸度闭合。 */
    [[nodiscard]] bool validatePhysicalFlashResult_(
        double pressure,
        double temperature,
        const Composition &z,
        const Result &result) const
    {
        if (!result.converged)
            return false;

        const double closureTolerance = std::max({
            10.0 * options_.fugacityTolerance,
            10.0 * options_.rachfordRiceTolerance,
            1.0e-10});

        double betaSum = 0.0;
        PhasePresence active(std::uint8_t{0});
        for (CompositionalPhase phase : phases_)
        {
            const std::size_t slot = static_cast<std::size_t>(phaseIndex(phase));
            const double beta = result.phaseMoleFraction[slot];
            if (!std::isfinite(beta) || beta < -closureTolerance)
                return false;
            if (beta > options_.phaseFractionTolerance)
            {
                active.add(phase);
                betaSum += beta;
                double xsum = 0.0;
                for (double x : result.composition[slot])
                {
                    if (!std::isfinite(x) || x < -closureTolerance)
                        return false;
                    xsum += x;
                }
                if (std::abs(xsum - 1.0) > closureTolerance)
                    return false;
            }
        }
        if (active.empty() || std::abs(betaSum - 1.0) > closureTolerance)
            return false;

        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            double reconstructed = 0.0;
            for (CompositionalPhase phase : phases_)
            {
                const std::size_t slot = static_cast<std::size_t>(phaseIndex(phase));
                if (active.contains(phase))
                    reconstructed += result.phaseMoleFraction[slot] *
                        result.composition[slot][c];
            }
            if (std::abs(reconstructed - z[c]) > closureTolerance)
                return false;
        }

        if (active.count() <= 1)
            return true;

        const CompositionalPhase reference = firstActive_(active);
        const std::size_t refSlot =
            static_cast<std::size_t>(phaseIndex(reference));
        const auto refThermo = phaseResultForThermodynamicRole_(
            reference, pressure, temperature, result.composition[refSlot]);
        for (CompositionalPhase phase : phases_)
        {
            if (!active.contains(phase) || phase == reference)
                continue;
            const std::size_t slot =
                static_cast<std::size_t>(phaseIndex(phase));
            const auto thermo = phaseResultForThermodynamicRole_(
                phase, pressure, temperature, result.composition[slot]);
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                // 组成被 active numerical floor 钳住时，该组分的等逸度方程已转成
                // 有界约束；继续比较 log(f_i) 只会把一个有效的 trace endpoint
                // 错判为失败。与已有 reference tests 一致，仅对非 floor-bound 组分终检。
                if (result.composition[refSlot][c] <=
                        10.0 * options_.compositionFloor ||
                    result.composition[slot][c] <=
                        10.0 * options_.compositionFloor)
                    continue;
                const double refActivity = std::max(
                    result.composition[refSlot][c] *
                        refThermo.fugacityCoefficient[c],
                    1.0e-300);
                const double activity = std::max(
                    result.composition[slot][c] *
                        thermo.fugacityCoefficient[c],
                    1.0e-300);
                if (std::abs(std::log(refActivity / activity)) >
                    10.0 * options_.fugacityTolerance)
                    return false;
            }
        }
        return true;
    }

    [[nodiscard]] Composition overallComposition_(const Result &result) const
    {
        Composition z{};
        for (std::size_t slot = 0; slot < 3; ++slot)
        {
            const double beta = result.phaseMoleFraction[slot];
            if (!(beta > options_.phaseFractionTolerance))
                continue;
            for (int i = 0; i < N; ++i)
                z[static_cast<std::size_t>(i)] +=
                    beta * result.composition[slot][static_cast<std::size_t>(i)];
        }
        normalize_(z);
        return z;
    }

    /** @brief 物理 O/G/W 角色调整后，在相数不变的前提下重新求解已选活动相集。 */
    [[nodiscard]] Result reequilibrateCanonicalRoles_(
        double pressure,
        double temperature,
        const Composition &z,
        Result seed) const
    {
        seed.presence = presenceFromFractions_(seed.phaseMoleFraction);
        const int count = seed.presence.count();
        if (count <= 1)
            return singlePhaseResult_(pressure, temperature, z, seed.presence);

        if (count == 2)
        {
            const auto pair = activePair_(seed.presence);
            const std::size_t first =
                static_cast<std::size_t>(phaseIndex(pair.first));
            const std::size_t second =
                static_cast<std::size_t>(phaseIndex(pair.second));
            Composition K{};
            for (int i = 0; i < N; ++i)
            {
                const std::size_t c = static_cast<std::size_t>(i);
                K[c] = std::clamp(
                    seed.composition[second][c] /
                        std::max(seed.composition[first][c], options_.compositionFloor),
                    1.0e-100, 1.0e100);
            }
            return flashTwo_(pressure, temperature, z, pair.first, pair.second, K);
        }

        deriveReferenceK_(seed);
        Result result = flashThree_(
            pressure, temperature, z, seed.vaporOilK, seed.waterOilK);
        if (result.converged)
            return result;

        return flashThree_(
            pressure, temperature, z,
            wilsonGasOilK_(pressure, temperature), initialWaterOilK_(z));
    }

    [[nodiscard]] bool canonicalizeNonAqueousRoles_(
        double pressure,
        double temperature,
        Result &result) const
    {
        constexpr std::size_t oil = 0;
        constexpr std::size_t gas = 1;
        const bool hasOil =
            result.phaseMoleFraction[oil] > options_.phaseFractionTolerance;
        const bool hasGas =
            result.phaseMoleFraction[gas] > options_.phaseFractionTolerance;
        if (!(hasOil && hasGas))
            return false;

        const Composition wilson = wilsonGasOilK_(pressure, temperature);
        double orientation = 0.0;
        double magnitude = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double lnK = std::log(std::clamp(wilson[c], 1.0e-30, 1.0e30));
            const double term =
                (result.composition[gas][c] - result.composition[oil][c]) * lnK;
            orientation += term;
            magnitude += std::abs(term);
        }

        constexpr double relativeOrientationTolerance = 1.0e-10;
        if (orientation >= -relativeOrientationTolerance * std::max(1.0, magnitude))
            return false;

        swapPhaseSlots_(result, oil, gas);
        return true;
    }

    [[nodiscard]] CompositionalPhase expectedSingleNonAqueousRole_(
        double pressure,
        double temperature,
        const Composition &composition) const
    {
        const Composition wilson = wilsonGasOilK_(pressure, temperature);
        double volatility = 0.0;
        double magnitude = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t c = static_cast<std::size_t>(i);
            const double term = composition[c] *
                std::log(std::clamp(wilson[c], 1.0e-30, 1.0e30));
            volatility += term;
            magnitude += std::abs(term);
        }

        constexpr double relativeVolatilityTolerance = 1.0e-10;
        return volatility >=
                -relativeVolatilityTolerance * std::max(1.0, magnitude)
            ? CompositionalPhase::Gas
            : CompositionalPhase::Oil;
    }

    [[nodiscard]] bool canonicalizeSingleNonAqueousRole_(
        double pressure,
        double temperature,
        Result &result) const
    {
        constexpr std::size_t oil = 0;
        constexpr std::size_t gas = 1;
        const bool hasOil =
            result.phaseMoleFraction[oil] > options_.phaseFractionTolerance;
        const bool hasGas =
            result.phaseMoleFraction[gas] > options_.phaseFractionTolerance;
        if (hasOil == hasGas)
            return false;

        const std::size_t active = hasOil ? oil : gas;
        const std::size_t expected = static_cast<std::size_t>(phaseIndex(
            expectedSingleNonAqueousRole_(
                pressure, temperature, result.composition[active])));
        if (active == expected)
            return false;

        swapPhaseSlots_(result, active, expected);
        result.presence = presenceFromFractions_(result.phaseMoleFraction);
        return true;
    }

    /** @brief SW 单一非水相在液/气根已合并时使用统一公共 O/G 标签。 */
    void canonicalizeCoalescedSwSingleNonAqueousRole_(
        double pressure,
        double temperature,
        Result &result) const
    {
        if (!eos_.usesSoreideWhitson())
            return;

        constexpr std::size_t oil = 0;
        constexpr std::size_t gas = 1;
        const bool hasOil =
            result.phaseMoleFraction[oil] > options_.phaseFractionTolerance;
        const bool hasGas =
            result.phaseMoleFraction[gas] > options_.phaseFractionTolerance;
        if (hasOil == hasGas)
            return;

        const std::size_t active = hasOil ? oil : gas;
        const auto &x = result.composition[active];
        const auto oilThermo = eos_.phaseResult(
            pressure, temperature, x, CompositionalPhase::Oil);
        const auto gasThermo = eos_.phaseResult(
            pressure, temperature, x, CompositionalPhase::Gas);
        const double rootScale = std::max({
            1.0,
            std::abs(oilThermo.compressibility),
            std::abs(gasThermo.compressibility)});
        constexpr double rootCoalescenceTolerance = 1.0e-8;
        if (std::abs(oilThermo.compressibility - gasThermo.compressibility) >
            rootCoalescenceTolerance * rootScale)
        {
            return;
        }

        (void)canonicalizeSingleNonAqueousRole_(
            pressure, temperature, result);
    }

    /**
     * @brief 折叠 SW 临界区内组成数值重合的重复活动相。
     *
     * 若两个活动相的全部组分都在显式容差内相同，它们不再代表可分辨的
     * 共存相。按相摩尔分率加权合并组成可严格保持总物料；W 参与时保留 W
     * 身份，O/G 重合时使用与单一非水相相同的 Wilson 挥发性规范。该判据
     * 比相分率截断严格得多，不会删除组成可分辨的真实小相。
     */
    void collapseCoincidentSwPhases_(
        double pressure,
        double temperature,
        Result &result) const
    {
        if (!eos_.usesSoreideWhitson() ||
            !(options_.coincidentPhaseCompositionTolerance > 0.0))
        {
            return;
        }

        constexpr std::array<std::array<std::size_t, 2>, 3> pairs{{
            {{0, 1}}, {{0, 2}}, {{1, 2}}
        }};
        bool changed = false;
        for (int pass = 0; pass < 2; ++pass)
        {
            bool mergedThisPass = false;
            for (const auto &pair : pairs)
            {
                const std::size_t a = pair[0];
                const std::size_t b = pair[1];
                const double betaA = result.phaseMoleFraction[a];
                const double betaB = result.phaseMoleFraction[b];
                if (!(betaA > options_.phaseFractionTolerance) ||
                    !(betaB > options_.phaseFractionTolerance))
                {
                    continue;
                }

                double maxDifference = 0.0;
                for (int i = 0; i < N; ++i)
                {
                    const std::size_t c = static_cast<std::size_t>(i);
                    maxDifference = std::max(
                        maxDifference,
                        std::abs(result.composition[a][c] -
                                 result.composition[b][c]));
                }
                if (maxDifference >
                    options_.coincidentPhaseCompositionTolerance)
                {
                    continue;
                }

                const double combinedBeta = betaA + betaB;
                Composition combined{};
                for (int i = 0; i < N; ++i)
                {
                    const std::size_t c = static_cast<std::size_t>(i);
                    combined[c] =
                        (betaA * result.composition[a][c] +
                         betaB * result.composition[b][c]) /
                        combinedBeta;
                }

                std::size_t destination = a;
                if (a == 2 || b == 2)
                {
                    destination = 2;
                }
                else
                {
                    destination = static_cast<std::size_t>(phaseIndex(
                        expectedSingleNonAqueousRole_(
                            pressure, temperature, combined)));
                }
                const std::size_t source = destination == a ? b : a;
                result.phaseMoleFraction[destination] = combinedBeta;
                result.phaseMoleFraction[source] = 0.0;
                result.composition[destination] = combined;
                result.composition[source] = combined;
                result.saturation[source] = 0.0;
                result.molarDensity[source] = 0.0;
                result.compressibility[source] = 1.0;
                mergedThisPass = true;
                changed = true;
                break;
            }
            if (!mergedThisPass)
                break;
        }

        if (!changed)
            return;
        result.presence = presenceFromFractions_(result.phaseMoleFraction);
        refreshCanonicalThermodynamics_(pressure, temperature, result);
    }

    [[nodiscard]] bool swThermodynamicRolesConsistent_(
        double pressure,
        double temperature,
        const Result &result,
        bool preserveRestrictedSingleNonAqueousRole) const
    {
        (void)pressure;
        if (!eos_.usesSoreideWhitson() || options_.waterComponent < 0)
            return true;

        const std::size_t wc = static_cast<std::size_t>(options_.waterComponent);
        int activeCount = 0;
        int aqueousSlot = -1;
        int oilSlot = -1;
        int gasSlot = -1;
        double richestWater = -1.0;
        int richestSlot = -1;

        for (int slot = 0; slot < 3; ++slot)
        {
            if (!(result.phaseMoleFraction[static_cast<std::size_t>(slot)] >
                  options_.phaseFractionTolerance))
                continue;
            ++activeCount;
            const auto role = thermodynamicRoles_[static_cast<std::size_t>(slot)];
            if (role == CompositionalPhase::Water) aqueousSlot = slot;
            if (role == CompositionalPhase::Oil) oilSlot = slot;
            if (role == CompositionalPhase::Gas) gasSlot = slot;
            const double xw = result.composition[static_cast<std::size_t>(slot)][wc];
            if (xw > richestWater)
            {
                richestWater = xw;
                richestSlot = slot;
            }
        }
        if (activeCount == 0)
            return false;

        // 物理：固定 Aqueous 角色必须对应活动相中最富 H2O 的相。
        // 不使用绝对含水阈值，避免残差在某个组成截断处发生离散模型切换。
        if (aqueousSlot >= 0)
        {
            constexpr double relTol = 1.0e-10;
            if (richestSlot != aqueousSlot &&
                richestWater > result.composition[static_cast<std::size_t>(aqueousSlot)][wc] +
                    relTol * std::max(1.0, std::abs(richestWater)))
                return false;
        }
        else if (richestSlot >= 0 &&
                 waterIsDominant_(result.composition[static_cast<std::size_t>(richestSlot)]))
        {
            return false;
        }

        // 物理：若只剩一个非水相，普通 unrestricted SW 路径继续使用历史的
        // 组成加权伪临界温度来规范公共 O/G 标签，保持既有 CO2-brine 与 P-T
        // 回归不变。restricted active-set 路径则额外尊重调用者明确指定的相集：
        // 若同一非水组成的 liquid/vapor cubic 根已经合并，EOS 本身已无法区分
        // 两个根，此时不应再用外部标签启发式否决一个已满足物料与逸度平衡的
        // G+W/O+W 解。两根仍明显分离时，restricted 路径也继续执行历史检查。
        if ((oilSlot >= 0) != (gasSlot >= 0))
        {
            const int slot = oilSlot >= 0 ? oilSlot : gasSlot;
            const auto &x = result.composition[static_cast<std::size_t>(slot)];

            bool rootsCoalesced = false;
            if (preserveRestrictedSingleNonAqueousRole)
            {
                const auto oilThermo = eos_.phaseResult(
                    pressure, temperature, x, CompositionalPhase::Oil);
                const auto gasThermo = eos_.phaseResult(
                    pressure, temperature, x, CompositionalPhase::Gas);
                const double rootScale = std::max({
                    1.0,
                    std::abs(oilThermo.compressibility),
                    std::abs(gasThermo.compressibility)});
                constexpr double rootCoalescenceTolerance = 1.0e-8;
                rootsCoalesced =
                    std::abs(oilThermo.compressibility - gasThermo.compressibility) <=
                    rootCoalescenceTolerance * rootScale;
            }

            if (!rootsCoalesced)
            {
                double pseudoTc = 0.0;
                for (int i = 0; i < N; ++i)
                    pseudoTc += x[static_cast<std::size_t>(i)] *
                        eos_.mixture().criticalTemperature(i);
                const double scale =
                    std::max({1.0, std::abs(temperature), std::abs(pseudoTc)});
                constexpr double relTol = 1.0e-8;
                if (std::abs(temperature - pseudoTc) > relTol * scale)
                {
                    const auto expected = temperature > pseudoTc
                        ? CompositionalPhase::Gas
                        : CompositionalPhase::Oil;
                    if (thermodynamicRoles_[static_cast<std::size_t>(slot)] != expected)
                        return false;
                }
            }
        }
        return true;
    }

    void canonicalizeSwThermodynamicRoles_(Result &result) const
    {
        Result canonical = result;
        canonical.phaseMoleFraction = {0.0, 0.0, 0.0};
        canonical.saturation = {0.0, 0.0, 0.0};
        canonical.compressibility = {1.0, 1.0, 1.0};
        canonical.molarDensity = {0.0, 0.0, 0.0};
        canonical.composition = {};

        for (int numericalSlot = 0; numericalSlot < 3; ++numericalSlot)
        {
            const auto role = thermodynamicRoles_[static_cast<std::size_t>(numericalSlot)];
            const std::size_t physicalSlot =
                static_cast<std::size_t>(phaseIndex(role));
            const std::size_t source = static_cast<std::size_t>(numericalSlot);
            canonical.phaseMoleFraction[physicalSlot] = result.phaseMoleFraction[source];
            canonical.saturation[physicalSlot] = result.saturation[source];
            canonical.composition[physicalSlot] = result.composition[source];
            canonical.compressibility[physicalSlot] = result.compressibility[source];
            canonical.molarDensity[physicalSlot] = result.molarDensity[source];
        }
        canonical.presence = presenceFromFractions_(canonical.phaseMoleFraction);
        deriveReferenceK_(canonical);
        result = canonical;
    }

    void canonicalizeFinal_(
        double pressure,
        double temperature,
        Result &result,
        bool preserveRestrictedSingleNonAqueousRole = false) const
    {
        if (!result.converged)
            return;

        if (eos_.usesSoreideWhitson())
        {
            if (!swThermodynamicRolesConsistent_(
                    pressure, temperature, result,
                    preserveRestrictedSingleNonAqueousRole))
            {
                result.converged = false;
                return;
            }
            canonicalizeSwThermodynamicRoles_(result);
            collapseCoincidentSwPhases_(pressure, temperature, result);
            // Restricted flashes carry an explicit active-set identity (O+W or G+W)
            // from the flow-state hysteresis logic.  Preserve that identity when the
            // cubic roots coalesce; the public Wilson label is only appropriate for
            // unrestricted equilibrium flashes.
            if (!preserveRestrictedSingleNonAqueousRole)
                canonicalizeCoalescedSwSingleNonAqueousRole_(
                    pressure, temperature, result);
            // Fixed SW thermodynamic roles determine the aqueous BIP matrix and
            // cubic root during the solve, but they do not by themselves define
            // the public Oil/Gas identity when both nonaqueous compositions are
            // in one-root regions.  Enforce the same volatility convention used
            // by PR/CPA after the role-to-public-slot mapping.  The complete
            // phase records move together, so no EOS role changes inside the
            // converged solve.
            (void)canonicalizeNonAqueousRoles_(
                pressure, temperature, result);
            result.presence = presenceFromFractions_(result.phaseMoleFraction);
            return;
        }

        // PR/CPA 保留既有全局 O/G/W 输出标准化；它们没有 SW 的“水相身份 -> BIP 矩阵”耦合，
        // 因而收敛后重排输出槽位不会改变 EOS 参数。
        constexpr int maximumRolePasses = 4;
        for (int pass = 0; pass < maximumRolePasses; ++pass)
        {
            result.presence = presenceFromFractions_(result.phaseMoleFraction);
            const Composition z = overallComposition_(result);
            const int waterSource = canonicalizeGlobalWaterRole_(result);

            const bool crossedGasWaterBoundary =
                waterSource == phaseIndex(CompositionalPhase::Gas);
            const bool crossedSingleNonAqueousBoundary =
                !preserveRestrictedSingleNonAqueousRole &&
                canonicalizeSingleNonAqueousRole_(
                    pressure, temperature, result);

            if (crossedGasWaterBoundary || crossedSingleNonAqueousBoundary)
            {
                result.presence = presenceFromFractions_(result.phaseMoleFraction);
                Result corrected = reequilibrateCanonicalRoles_(
                    pressure, temperature, z, result);
                if (!corrected.converged)
                {
                    result.converged = false;
                    return;
                }
                result = corrected;
                continue;
            }

            (void)canonicalizeNonAqueousRoles_(
                pressure, temperature, result);
            result.presence = presenceFromFractions_(result.phaseMoleFraction);
            refreshCanonicalThermodynamics_(pressure, temperature, result);
            return;
        }

        result.converged = false;
    }

    [[nodiscard]] PhasePresence presenceFromFractions_(const std::array<double, 3> &beta) const
    {
        PhasePresence presence(0);
        for (int p = 0; p < 3; ++p)
        {
            if (beta[static_cast<std::size_t>(p)] > options_.phaseFractionTolerance)
                presence.add(phases_[static_cast<std::size_t>(p)]);
        }
        return presence;
    }

    [[nodiscard]] static CompositionalPhase firstActive_(PhasePresence active)
    {
        for (CompositionalPhase p : phases_)
            if (active.contains(p)) return p;
        throw std::logic_error("PhasePresence is empty.");
    }

    [[nodiscard]] static std::pair<CompositionalPhase, CompositionalPhase>
    activePair_(PhasePresence active)
    {
        std::array<CompositionalPhase, 2> pair{};
        int n = 0;
        for (CompositionalPhase p : phases_)
            if (active.contains(p) && n < 2) pair[static_cast<std::size_t>(n++)] = p;
        if (n != 2) throw std::logic_error("activePair requires exactly two phases.");
        return {pair[0], pair[1]};
    }

