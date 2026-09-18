// Pointwise Figure-7 WLV-WL verification.
// Reuse the exact frozen Jia/Amani benchmark factory and boundary routines.
#define main athabasca_reference_entrypoint
#include "main.cpp"
#undef main

namespace
{
constexpr double sideOffsetMPa = 0.02;
constexpr double boundaryPathToleranceMPa = 1.0e-4;
constexpr double compositionPathTolerance = 1.0e-6;
constexpr double rootTolerance = 1.0e-8;
constexpr double distinctRootTolerance = 1.0e-6;
constexpr double gibbsTolerance = 1.0e-8;

const auto oilWaterPresence()
{
    return MPMC::PhasePresence(static_cast<std::uint8_t>(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit));
}

Flash::Result continuedOwAt(
    const Flash &flash,
    double temperatureK,
    const Composition &z,
    double targetPressureMPa)
{
    const auto ow = oilWaterPresence();
    double p = 30.0;
    auto state = flash.flashRestricted(p * 1.0e6, temperatureK, z, ow);
    if (!state.converged)
        return state;
    while (p > targetPressureMPa)
    {
        p = std::max(targetPressureMPa, p - 0.05);
        state = flash.flashRestricted(
            p * 1.0e6, temperatureK, z, ow, state.composition);
        if (!state.converged)
            return state;
    }
    return state;
}

double compositionGap(
    const Flash::Result &a,
    const Flash::Result &b)
{
    if (!a.converged || !b.converged ||
        a.presence.bits() != b.presence.bits())
        return std::numeric_limits<double>::infinity();

    double gap = 0.0;
    for (std::size_t phase = 0; phase < 3; ++phase)
    {
        const auto role = static_cast<MPMC::CompositionalPhase>(phase);
        if (!a.presence.contains(role))
            continue;
        for (std::size_t i = 0; i < a.composition[phase].size(); ++i)
            gap = std::max(
                gap,
                std::abs(a.composition[phase][i] - b.composition[phase][i]));
    }
    return gap;
}

double compositionL1(const Composition &a, const Composition &b)
{
    double value = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        value += std::abs(a[i] - b[i]);
    return value;
}

double dimensionlessGibbs(
    const Eos &eos,
    double pressureMPa,
    double temperatureK,
    const Flash::Result &result)
{
    if (!result.converged)
        return std::numeric_limits<double>::quiet_NaN();

    const double pressure = pressureMPa * 1.0e6;
    double g = 0.0;
    for (std::size_t phase = 0; phase < 3; ++phase)
    {
        const auto role = static_cast<MPMC::CompositionalPhase>(phase);
        if (!result.presence.contains(role))
            continue;
        const double beta = result.phaseMoleFraction[phase];
        if (!(beta > 0.0))
            continue;

        // Recompute on the active role's fixed root: Oil/Water liquid,
        // Gas vapor. Candidate-root minimization is not allowed to rewrite
        // an already active phase.
        const auto thermo = eos.phaseResult(
            pressure, temperatureK, result.composition[phase], role, false);
        for (std::size_t i = 0; i < result.composition[phase].size(); ++i)
        {
            const double x = result.composition[phase][i];
            if (!(x > 0.0))
                continue;
            const double phi = thermo.fugacityCoefficient[i];
            if (!(phi > 0.0) || !std::isfinite(phi))
                return std::numeric_limits<double>::quiet_NaN();
            g += beta * x * (std::log(x) + std::log(phi));
        }
    }
    return g;
}

bool certified(
    const Flash &flash,
    double pressureMPa,
    double temperatureK,
    const Composition &z,
    const Flash::Result &result)
{
    if (!result.converged || maxMaterialClosure(z, result) > 1.0e-8)
        return false;
    const auto stability = flash.stabilityTest(
        pressureMPa * 1.0e6, temperatureK, z,
        result.presence, result.composition);
    return stability.valid && stability.stable;
}

struct SideCertificate
{
    bool owConverged{false};
    bool stabilityValid{false};
    bool gasUnstable{false};
    bool onlyGasMissingDirection{false};
    double gasTrialSum{std::numeric_limits<double>::quiet_NaN()};
    Composition incipientGas{};
    bool incipientFinite{false};
    double selectedZ{std::numeric_limits<double>::quiet_NaN()};
    double vaporZ{std::numeric_limits<double>::quiet_NaN()};
    double liquidZ{std::numeric_limits<double>::quiet_NaN()};
    bool selectedIsVapor{false};
    bool vaporLiquidDistinct{false};
};

SideCertificate inspectRestrictedOw(
    const Eos &eos,
    const Flash &flash,
    double pressureMPa,
    double temperatureK,
    const Composition &z,
    const Flash::Result &ow)
{
    SideCertificate c;
    const auto expected = oilWaterPresence();
    c.owConverged = ow.converged &&
        ow.presence.bits() == expected.bits() &&
        maxMaterialClosure(z, ow) <= 1.0e-8;
    if (!c.owConverged)
        return c;

    const auto stability = flash.stabilityTest(
        pressureMPa * 1.0e6, temperatureK, z,
        ow.presence, ow.composition);
    c.stabilityValid = stability.valid;
    if (!stability.valid)
        return c;

    const std::size_t gas = static_cast<std::size_t>(
        MPMC::phaseIndex(MPMC::CompositionalPhase::Gas));
    const std::size_t oil = static_cast<std::size_t>(
        MPMC::phaseIndex(MPMC::CompositionalPhase::Oil));
    const std::size_t waterPhase = static_cast<std::size_t>(
        MPMC::phaseIndex(MPMC::CompositionalPhase::Water));

    c.gasUnstable = stability.missingPhaseUnstable[gas];
    c.onlyGasMissingDirection =
        !stability.missingPhaseUnstable[oil] &&
        !stability.missingPhaseUnstable[waterPhase];
    c.gasTrialSum = stability.trialSum[gas];
    c.incipientGas = stability.incipientComposition[gas];

    double sum = 0.0;
    c.incipientFinite = true;
    for (double x : c.incipientGas)
    {
        c.incipientFinite =
            c.incipientFinite && std::isfinite(x) && x >= 0.0;
        sum += x;
    }
    c.incipientFinite =
        c.incipientFinite && std::abs(sum - 1.0) <= 1.0e-8;
    if (!c.incipientFinite || !c.gasUnstable)
        return c;

    try
    {
        const double pressure = pressureMPa * 1.0e6;
        // This is the missing-phase candidate policy used by the Jia audit.
        const auto selected = eos.phaseResult(
            pressure, temperatureK, c.incipientGas,
            MPMC::CompositionalPhase::Gas, true);
        const auto vapor = eos.phaseResult(
            pressure, temperatureK, c.incipientGas,
            MPMC::CompositionalPhase::Gas, false);
        const auto liquid = eos.phaseResult(
            pressure, temperatureK, c.incipientGas,
            MPMC::CompositionalPhase::Oil, false);
        c.selectedZ = selected.compressibility;
        c.vaporZ = vapor.compressibility;
        c.liquidZ = liquid.compressibility;
        const double scale = std::max(1.0, std::abs(c.vaporZ));
        c.selectedIsVapor =
            std::abs(c.selectedZ - c.vaporZ) <= rootTolerance * scale;
        c.vaporLiquidDistinct =
            std::abs(c.vaporZ - c.liquidZ) > distinctRootTolerance;
    }
    catch (const std::exception &)
    {
        c.selectedIsVapor = false;
        c.vaporLiquidDistinct = false;
    }
    return c;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 2)
            throw std::invalid_argument(
                "usage: figure7_pointwise_probe OUTPUT_DIR");

        const std::filesystem::path out = argv[1];
        std::filesystem::create_directories(out);
        std::ofstream rows(out / "figure7_pointwise.csv");
        std::ofstream gates(out / "figure7_pointwise_gate.csv");
        if (!rows || !gates)
            throw std::runtime_error("Cannot open Figure-7 pointwise outputs.");

        rows << std::scientific << std::setprecision(12);
        gates << std::scientific << std::setprecision(12);
        rows
            << "T_K,experimental_P_MPa,independent_P_MPa,continuation_P_MPa,"
               "path_gap_MPa,independent_transition_count,"
               "continuation_transition_count,probe_low_P_MPa,probe_high_P_MPa,"
               "cold_warm_gap_low,cold_warm_gap_high,"
               "low_ow_gas_unstable,low_only_gas_missing_direction,"
               "low_gas_trial_sum,low_incipient_xH2O,"
               "low_candidate_selected_Z,low_candidate_vapor_Z,"
               "low_candidate_liquid_Z,low_candidate_selected_is_vapor,"
               "low_candidate_vapor_liquid_distinct,"
               "low_global_phase_code,low_global_beta_gas,"
               "low_global_certified,low_global_gas_Z,"
               "low_global_gas_vapor_Z,low_global_gas_root_match,"
               "low_incipient_global_gas_L1,low_global_gibbs_RT,"
               "low_restricted_ow_gibbs_RT,"
               "high_ow_gas_unstable,high_global_phase_code,"
               "high_global_beta_gas,high_global_certified,"
               "path_pass,branch_pass,topology_pass,gas_identity_pass,"
               "gibbs_pass,point_pass\n";
        gates
            << "T_K,status,path_gap_MPa,cold_warm_gap_max,"
               "low_phase_code,high_phase_code,gas_selected_is_vapor,"
               "global_minus_restricted_gibbs_RT,criterion\n";

        Eos eos = makeJiaCase1Cpa();
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = static_cast<int>(water);
        options.maximumIterations = 240;
        options.maximumStabilityIterations = 160;
        options.cpaSelectGibbsMinimumRoot = true;
        Flash flash(eos, options);

        // Authoritative Amani source: 55.9 wt% bitumen + 44.1 wt% water.
        const Composition z = feedFromWaterMassFraction(0.441);
        const std::array<std::array<double, 2>, 6> sourcePoints{{
            {{522.9, 4.2}},
            {{573.3, 9.0}},
            {{583.0, 10.8}},
            {{593.0, 12.8}},
            {{603.6, 16.2}},
            {{613.2, 20.1}}
        }};

        bool allPass = true;
        for (const auto &point : sourcePoints)
        {
            const double t = point[0];
            const double pExp = point[1];
            const auto independent = findWlvWlBoundary(flash, t, z);
            const auto continued =
                findWlvWlBoundaryContinuation(flash, t, z);

            const double pathGap =
                independent.found && continued.found
                ? std::abs(independent.pressureMPa - continued.pressureMPa)
                : std::numeric_limits<double>::infinity();
            const bool pathPass =
                independent.found && continued.found &&
                independent.transitionCount == 1 &&
                continued.transitionCount == 1 &&
                pathGap <= boundaryPathToleranceMPa;

            const double boundaryP =
                independent.found && continued.found
                ? 0.5 * (independent.pressureMPa + continued.pressureMPa)
                : std::numeric_limits<double>::quiet_NaN();
            const double pLow = boundaryP - sideOffsetMPa;
            const double pHigh = boundaryP + sideOffsetMPa;

            const auto ow = oilWaterPresence();
            const auto coldLow = flash.flashRestricted(
                pLow * 1.0e6, t, z, ow);
            const auto warmLow = continuedOwAt(flash, t, z, pLow);
            const auto coldHigh = flash.flashRestricted(
                pHigh * 1.0e6, t, z, ow);
            const auto warmHigh = continuedOwAt(flash, t, z, pHigh);
            const double lowGap = compositionGap(coldLow, warmLow);
            const double highGap = compositionGap(coldHigh, warmHigh);
            const bool branchPass =
                std::isfinite(lowGap) && std::isfinite(highGap) &&
                lowGap <= compositionPathTolerance &&
                highGap <= compositionPathTolerance;

            const auto lowOw = inspectRestrictedOw(
                eos, flash, pLow, t, z, warmLow);
            const auto highOw = inspectRestrictedOw(
                eos, flash, pHigh, t, z, warmHigh);

            const auto globalLow = flash.flash(pLow * 1.0e6, t, z);
            const auto globalHigh = flash.flash(pHigh * 1.0e6, t, z);
            const bool lowCertified =
                certified(flash, pLow, t, z, globalLow);
            const bool highCertified =
                certified(flash, pHigh, t, z, globalHigh);

            const auto all = MPMC::PhasePresence::all();
            const bool lowThreePhase =
                globalLow.converged &&
                globalLow.presence.bits() == all.bits();
            const bool highTwoLiquid =
                globalHigh.converged &&
                globalHigh.presence.bits() == ow.bits();
            const std::size_t gas = static_cast<std::size_t>(
                MPMC::phaseIndex(MPMC::CompositionalPhase::Gas));
            const double lowBetaGas = lowThreePhase
                ? globalLow.phaseMoleFraction[gas] : 0.0;
            const double highBetaGas =
                globalHigh.converged ? globalHigh.phaseMoleFraction[gas] : 0.0;

            double lowGlobalGasZ =
                std::numeric_limits<double>::quiet_NaN();
            double lowGlobalVaporZ =
                std::numeric_limits<double>::quiet_NaN();
            bool lowGlobalGasRootMatch = false;
            double incipientGlobalGasL1 =
                std::numeric_limits<double>::quiet_NaN();
            if (lowThreePhase)
            {
                const auto vapor = eos.phaseResult(
                    pLow * 1.0e6, t, globalLow.composition[gas],
                    MPMC::CompositionalPhase::Gas, false);
                lowGlobalGasZ = globalLow.compressibility[gas];
                lowGlobalVaporZ = vapor.compressibility;
                lowGlobalGasRootMatch =
                    std::abs(lowGlobalGasZ - lowGlobalVaporZ) <=
                    rootTolerance * std::max(1.0, std::abs(lowGlobalVaporZ));
                if (lowOw.incipientFinite)
                    incipientGlobalGasL1 = compositionL1(
                        lowOw.incipientGas, globalLow.composition[gas]);
            }

            const bool topologyPass =
                lowOw.owConverged && highOw.owConverged &&
                lowOw.stabilityValid && highOw.stabilityValid &&
                lowOw.gasUnstable && !highOw.gasUnstable &&
                lowOw.onlyGasMissingDirection &&
                lowThreePhase && highTwoLiquid &&
                lowCertified && highCertified &&
                lowBetaGas > options.phaseFractionTolerance &&
                highBetaGas <= options.phaseFractionTolerance;

            // A vapor-like incipient phase need not have two distinct density
            // roots. In a one-root region the Gas-role vapor evaluation and
            // liquid evaluation legitimately coincide. Phase identity is
            // certified by the nontrivial missing-Gas TPD direction plus the
            // released global Gas phase retaining its Gas-role root.
            const bool gasIdentityPass =
                lowOw.incipientFinite &&
                lowOw.gasUnstable &&
                lowOw.onlyGasMissingDirection &&
                lowOw.selectedIsVapor &&
                lowGlobalGasRootMatch;

            const double gGlobalLow =
                dimensionlessGibbs(eos, pLow, t, globalLow);
            const double gOwLow =
                dimensionlessGibbs(eos, pLow, t, warmLow);
            const bool gibbsPass =
                std::isfinite(gGlobalLow) && std::isfinite(gOwLow) &&
                gGlobalLow <= gOwLow + gibbsTolerance;

            const bool pass =
                pathPass && branchPass && topologyPass &&
                gasIdentityPass && gibbsPass;
            allPass = allPass && pass;

            rows << t << ',' << pExp << ','
                 << independent.pressureMPa << ','
                 << continued.pressureMPa << ',' << pathGap << ','
                 << independent.transitionCount << ','
                 << continued.transitionCount << ','
                 << pLow << ',' << pHigh << ','
                 << lowGap << ',' << highGap << ','
                 << lowOw.gasUnstable << ','
                 << lowOw.onlyGasMissingDirection << ','
                 << lowOw.gasTrialSum << ','
                 << (lowOw.incipientFinite ? lowOw.incipientGas[water] :
                     std::numeric_limits<double>::quiet_NaN()) << ','
                 << lowOw.selectedZ << ',' << lowOw.vaporZ << ','
                 << lowOw.liquidZ << ',' << lowOw.selectedIsVapor << ','
                 << lowOw.vaporLiquidDistinct << ','
                 << (globalLow.converged ? globalLow.presence.bits() : 0)
                 << ',' << lowBetaGas << ',' << lowCertified << ','
                 << lowGlobalGasZ << ',' << lowGlobalVaporZ << ','
                 << lowGlobalGasRootMatch << ','
                 << incipientGlobalGasL1 << ','
                 << gGlobalLow << ',' << gOwLow << ','
                 << highOw.gasUnstable << ','
                 << (globalHigh.converged ? globalHigh.presence.bits() : 0)
                 << ',' << highBetaGas << ',' << highCertified << ','
                 << pathPass << ',' << branchPass << ','
                 << topologyPass << ',' << gasIdentityPass << ','
                 << gibbsPass << ',' << pass << '\n';

            gates << t << ',' << (pass ? "PASS" : "FAIL") << ','
                  << pathGap << ',' << std::max(lowGap, highGap) << ','
                  << (globalLow.converged ? globalLow.presence.bits() : 0)
                  << ','
                  << (globalHigh.converged ? globalHigh.presence.bits() : 0)
                  << ',' << lowOw.selectedIsVapor << ','
                  << (gGlobalLow - gOwLow)
                  << ",independent_and_continued_boundary_agree;"
                     "low_side_is_certified_OGW;high_side_is_certified_OW;"
                     "missing_direction_is_Gas_and_candidate_uses_Gas_vapor_or_unique_root;"
                     "released_state_has_no_higher_G\n";

            std::cout << "Figure7 point T=" << t
                      << " K path_gap=" << pathGap
                      << " MPa low_code="
                      << (globalLow.converged ? globalLow.presence.bits() : 0)
                      << " high_code="
                      << (globalHigh.converged ? globalHigh.presence.bits() : 0)
                      << " gas_vapor=" << lowOw.selectedIsVapor
                      << " status=" << (pass ? "PASS" : "FAIL")
                      << std::endl;
        }

        std::cout << (allPass
            ? "FIGURE7_POINTWISE_BOUNDARY_PASS"
            : "FIGURE7_POINTWISE_BOUNDARY_BLOCKED") << std::endl;
        return allPass ? 0 : 2;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Figure7 pointwise probe failed: "
                  << e.what() << std::endl;
        return 1;
    }
}
