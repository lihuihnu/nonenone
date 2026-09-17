#!/usr/bin/env python3
from pathlib import Path


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one anchor, found {count}")
    path.write_text(text.replace(old, new, 1))


root = Path(__file__).resolve().parents[1]
flash = root / "models/include/natural/thermo/three_phase_flash.hpp"
test = root / "test/src/unit/sw_flash_recovery_test.cpp"

old_selection = r'''            const int p = phaseIndex(candidate);
            Trial best;
            best.sum = -std::numeric_limits<double>::infinity();

            const auto seeds = stabilitySeeds_(
                candidate, pressure, temperature, z, referenceComposition);
            for (std::size_t seedIndex = 0; seedIndex < seeds.size; ++seedIndex)
            {
                const auto &seed = seeds.values[seedIndex];
                const Trial trial = stabilityTrial_(
                    candidate, pressure, temperature, logReference,
                    referenceComposition, seed);
                if (trial.valid && (!best.valid || trial.sum > best.sum))
                    best = trial;
            }

            if (!best.valid)
            {
                result.valid = false;
                result.stable = false; // 数值：试探失败时按“不稳定”保守处理，避免错误抑制新相出现。
                continue;
            }
'''

new_selection = r'''            const int p = phaseIndex(candidate);
            Trial best;
            best.sum = -std::numeric_limits<double>::infinity();
            Trial bestWrongSwRole;
            bestWrongSwRole.sum = -std::numeric_limits<double>::infinity();

            const auto seeds = stabilitySeeds_(
                candidate, pressure, temperature, z, referenceComposition);
            for (std::size_t seedIndex = 0; seedIndex < seeds.size; ++seedIndex)
            {
                const auto &seed = seeds.values[seedIndex];
                const Trial trial = stabilityTrial_(
                    candidate, pressure, temperature, logReference,
                    referenceComposition, seed);
                if (!trial.valid)
                    continue;

                // SW fixes the thermodynamic role (including the aqueous BIP
                // matrix) for the duration of each TPD solve. A non-aqueous
                // Oil/Gas trial can nevertheless converge to a strongly
                // H2O-rich stationary composition that the public phase model
                // assigns to Water. That stationary point was evaluated with
                // the wrong SW role and therefore cannot certify appearance of
                // a non-aqueous phase. Ignore it for this candidate and let the
                // independently solved Water candidate decide aqueous
                // stability with the proper aqueous BIP model. Filter before
                // selecting the strongest multi-start stationary point so a
                // wrong-role basin cannot mask a weaker, genuine Gas/Oil one.
                const bool wrongSwRole =
                    eos_.usesSoreideWhitson() &&
                    candidate != CompositionalPhase::Water &&
                    waterIsDominant_(trial.composition);
                if (wrongSwRole)
                {
                    if (!bestWrongSwRole.valid ||
                        trial.sum > bestWrongSwRole.sum)
                    {
                        bestWrongSwRole = trial;
                    }
                    continue;
                }

                if (!best.valid || trial.sum > best.sum)
                    best = trial;
            }

            if (!best.valid)
            {
                if (bestWrongSwRole.valid)
                {
                    // A stationary solve did converge, but only in a physical
                    // role excluded from this candidate. Preserve it for
                    // diagnostics without turning a role-classification issue
                    // into an invalid stability calculation.
                    result.trialSum[static_cast<std::size_t>(p)] =
                        bestWrongSwRole.sum;
                    result.incipientComposition[static_cast<std::size_t>(p)] =
                        bestWrongSwRole.composition;
                    result.missingPhaseUnstable[static_cast<std::size_t>(p)] = false;
                    continue;
                }

                result.valid = false;
                result.stable = false; // 数值：试探失败时按“不稳定”保守处理，避免错误抑制新相出现。
                continue;
            }
'''
replace_once(flash, old_selection, new_selection)

insert_anchor = r'''
} // namespace

int main()
'''
new_regression = r'''
Eos makeLmhBoundarySw()
{
    constexpr std::array<double, 4> lmhTc{
        647.30, 354.1916431226766, 605.78, 751.00};
    constexpr std::array<double, 4> lmhPc{
        22.048e6, 4.065799256505576e6, 2.175e6, 1.654e6};
    constexpr std::array<double, 4> lmhVc{
        5.594803743e-5, 1.95564637471291e-4,
        6.252498825299838e-4, 1.0193008374141067e-3};
    constexpr std::array<double, 4> lmhOmega{
        0.344, 0.1498936802973978, 0.618, 0.957};
    constexpr std::array<double, 4> lmhMw{
        0.018015, 0.04606110037174722, 0.140960, 0.280990};
    constexpr std::array<std::array<double, 4>, 4> lmhKij{{
        {{0.0, 0.5, 0.5, 0.5}},
        {{0.5, 0.0, 0.0, 0.0}},
        {{0.5, 0.0, 0.0, 0.0}},
        {{0.5, 0.0, 0.0, 0.0}}
    }};

    MPMC::CompositionalMixture<Indices> mixture(
        lmhTc, lmhPc, lmhVc, lmhOmega, lmhMw, lmhKij);
    Eos eos(
        0.45724, 0.07780, std::move(mixture), 1,
        2.4142135623730951, -0.4142135623730951, 1.0e-30);

    Eos::SoreideWhitsonOptions sw;
    sw.waterComponent = 0;
    sw.salinityMolality = 0.0;
    sw.aqueousWaterBip[0] = [](double, double) { return 0.0; };
    for (int component = 1; component < 4; ++component)
    {
        const std::size_t c = static_cast<std::size_t>(component);
        const double componentTc = lmhTc[c];
        const double componentOmega = lmhOmega[c];
        sw.aqueousWaterBip[c] =
            [componentTc, componentOmega](double temperature, double salinity) {
                return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
                    temperature, componentTc, componentOmega, salinity);
            };
    }
    eos.configureSoreideWhitson(std::move(sw));
    eos.configureAqueousCompositionDomain(0, 0.02);
    return eos;
}

void checkLmhSinglePhaseRejectsWrongRoleGasTpdBasin()
{
    // Real cell-540 boundary state from the 60x20 LMH SW displacement. The
    // pre-fix Gas multi-start search selected a 99.48 mol% H2O stationary point
    // (trialSum > 1) even though that composition belongs to the explicit
    // aqueous physical-property domain. The independently solved Water trial
    // lies outside that domain, so the physical state is still Oil-only.
    constexpr double pressure = 278.772e5;
    constexpr double temperature = 653.2;
    constexpr Composition z{
        0.91456, 0.0385195, 0.0315411, 0.0153792};

    const auto eos = makeLmhBoundarySw();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    const Flash flash(eos, options);
    const std::array<Composition, 3> compositions{z, z, z};
    const auto stability = flash.stabilityTest(
        pressure,
        temperature,
        z,
        MPMC::PhasePresence::oilOnly(),
        compositions);

    require(stability.valid,
            "LMH SW Oil-only boundary stability must remain numerically valid");
    require(stability.stable,
            "LMH SW Oil-only boundary must not create a wrong-role phase");
    require(!stability.missingPhaseUnstable[1],
            "water-like nonaqueous TPD basin must not appear as Gas");
    require(!stability.missingPhaseUnstable[2],
            "unsupported aqueous TPD basin must not appear as Water");
    require(!eos.aqueousVolumeCompositionSupported(
                stability.incipientComposition[1]),
            "selected Gas TPD stationary point must remain outside Water role domain");
}

'''
replace_once(test, insert_anchor, new_regression + insert_anchor)
replace_once(
    test,
    "        checkPressureContinuationThreePhase(eos);\n",
    "        checkPressureContinuationThreePhase(eos);\n"
    "        checkLmhSinglePhaseRejectsWrongRoleGasTpdBasin();\n")

print("SW role-aware TPD candidate patch applied")
