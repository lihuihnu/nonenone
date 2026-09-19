// Focused verification, not parameter regression. Reuse the exact benchmark
// factory in one translation unit so no Athabasca or OIL_HEAVY values drift.
#define main athabasca_reference_entrypoint
#include "main.cpp"
#undef main

namespace {
struct ProbeCheck {
    bool valid{false};
    bool stable{false};
    double closure{0.0};
    double fugacityNorm{0.0};
    double rootMismatch{0.0};
    double gibbsRT{0.0};
};

ProbeCheck recordProbe(const Eos &eos, const Flash &flash,
    double t, double pMPa, double w, const Composition &z,
    const std::string &method, const Flash::Result &r,
    std::ostream &states, std::ostream &phases)
{
    const double p = pMPa * 1.0e6;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    ProbeCheck c;
    c.closure = c.fugacityNorm = c.rootMismatch = c.gibbsRT = nan;
    if (r.converged) {
        c.closure = maxMaterialClosure(z, r);
        const auto st = flash.stabilityTest(p, t, z, r.presence, r.composition);
        c.valid = st.valid;
        c.stable = st.stable;
        c.fugacityNorm = c.rootMismatch = c.gibbsRT = 0.0;
        std::array<double, 5> reference{};
        bool haveReference = false;
        for (std::size_t k = 0; k < 3; ++k) {
            const auto role = static_cast<MPMC::CompositionalPhase>(k);
            if (!r.presence.contains(role)) continue;
            const auto &x = r.composition[k];
            // Recompute on the stored role's fixed root, independently of
            // the candidate-root Gibbs option used by stability searches.
            const auto selected = eos.phaseResult(p, t, x, role, false);
            const auto liquid = eos.phaseResult(p, t, x, true, false);
            const auto vapor = eos.phaseResult(p, t, x, false, false);
            c.rootMismatch = std::max(c.rootMismatch,
                std::abs(r.compressibility[k] - selected.compressibility));
            for (std::size_t i = 0; i < x.size(); ++i) {
                const double f = selected.fugacity[i];
                if (!(f > 0.0) || !std::isfinite(f)) {
                    c.valid = false;
                    c.fugacityNorm = std::numeric_limits<double>::infinity();
                    continue;
                }
                const double logf = std::log(f);
                if (haveReference) c.fugacityNorm = std::max(
                    c.fugacityNorm, std::abs(logf - reference[i]));
                else reference[i] = logf;
                if (x[i] > 0.0) c.gibbsRT += r.phaseMoleFraction[k] *
                    x[i] * std::log(x[i] * selected.fugacityCoefficient[i]);
            }
            haveReference = true;
            phases << t << ',' << pMPa << ',' << w << ',' << method << ','
                << k << ',' << r.phaseMoleFraction[k] << ',' << x[0] << ','
                << r.compressibility[k] << ',' << selected.compressibility
                << ',' << liquid.compressibility << ',' << vapor.compressibility
                << ',' << p / (MPMC::units::gasConstant * t * selected.compressibility)
                << '\n';
        }
        c.valid = c.valid && std::isfinite(c.closure) && c.closure <= 1.0e-8 &&
            std::isfinite(c.fugacityNorm) && c.fugacityNorm <= 1.0e-6 &&
            std::isfinite(c.rootMismatch) && c.rootMismatch <= 1.0e-8;
    }
    states << t << ',' << pMPa << ',' << w << ',' << method << ','
        << r.converged << ',' << (r.converged ? r.presence.bits() : 0) << ','
        << c.valid << ',' << c.stable << ',' << c.closure << ','
        << c.fugacityNorm << ',' << c.rootMismatch << ',' << c.gibbsRT << '\n';
    return c;
}

double compositionGap(const Flash::Result &a, const Flash::Result &b) {
    if (!a.converged || !b.converged || a.presence.bits() != b.presence.bits())
        return std::numeric_limits<double>::infinity();
    double gap = 0.0;
    for (std::size_t k = 0; k < 3; ++k) {
        if (!a.presence.contains(static_cast<MPMC::CompositionalPhase>(k))) continue;
        for (std::size_t i = 0; i < 5; ++i)
            gap = std::max(gap, std::abs(a.composition[k][i] - b.composition[k][i]));
    }
    return gap;
}
}

int main(int argc, char **argv) {
    try {
        if (argc != 2) throw std::invalid_argument("usage: branch_probe OUTPUT_DIR");
        const std::filesystem::path out = argv[1];
        std::filesystem::create_directories(out);
        std::ofstream states(out / "state_summary.csv");
        std::ofstream phases(out / "phase_probe.csv");
        std::ofstream gates(out / "path_gate.csv");
        if (!states || !phases || !gates) throw std::runtime_error("Cannot open probe output");
        states << std::setprecision(17);
        phases << std::setprecision(17);
        gates << std::setprecision(17);
        states << "T_K,P_MPa,water_mass_fraction,method,converged,phase_code,certificate_valid,globally_stable,material_closure,selected_root_log_fugacity_spread,stored_selected_Z_mismatch,gibbs_RT\n";
        phases << "T_K,P_MPa,water_mass_fraction,method,phase_slot,beta,xH2O,Z_stored,Z_selected,Z_liquid,Z_vapor,selected_molar_density\n";
        gates << "T_K,P_MPa,water_mass_fraction,status,cold_warm_composition_gap,global_warm_composition_gap,criterion\n";
        Eos eos = makeJiaCase1Cpa();
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = 0;
        options.maximumIterations = 240;
        options.maximumStabilityIterations = 160;
        options.cpaSelectGibbsMinimumRoot = true;
        Flash flash(eos, options);
        const auto ow = MPMC::PhasePresence(static_cast<std::uint8_t>(
            MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit));
        const std::array<std::array<double, 3>, 5> probes{{
            {{548.2, 6.91, 0.559}}, {{573.1, 9.52, 0.559}},
            {{573.1, 9.52, 0.441}}, {{583.0, 11.0, 0.441}},
            {{603.5, 15.32, 0.559}}
        }};
        bool allPass = true;
        for (const auto &point : probes) {
            const double t = point[0], p = point[1], w = point[2];
            const Composition z = feedFromWaterMassFraction(w);
            const auto global = flash.flash(p * 1.0e6, t, z);
            const auto cold = flash.flashRestricted(p * 1.0e6, t, z, ow);
            double walkP = 30.0;
            auto warm = flash.flashRestricted(walkP * 1.0e6, t, z, ow);
            while (warm.converged && walkP > p) {
                walkP = std::max(p, walkP - 0.5);
                warm = flash.flashRestricted(walkP * 1.0e6, t, z, ow, warm.composition);
            }
            const auto cg = recordProbe(eos, flash, t, p, w, z, "UNRESTRICTED", global, states, phases);
            const auto cc = recordProbe(eos, flash, t, p, w, z, "COLD_OW", cold, states, phases);
            const auto cw = recordProbe(eos, flash, t, p, w, z, "CONTINUED_OW", warm, states, phases);
            const bool owGlobal = global.converged && global.presence.bits() == ow.bits();
            const double gap = compositionGap(cold, warm);
            const double globalGap = compositionGap(global, warm);
            // A metastable OW branch is legitimate below a WLV boundary, but
            // it cannot validate OW equilibrium. Only require all paths to
            // agree when unrestricted equilibrium actually has two liquids.
            // Even a metastable liquid branch must remain root-consistent,
            // fugacity-closed and seed-independent. A global three-phase
            // equilibrium must have no higher G than the restricted branch.
            const bool pass = cg.valid && cg.stable && cc.valid && cw.valid &&
                gap <= 1.0e-6 && cg.gibbsRT <= cc.gibbsRT + 1.0e-8 &&
                cg.gibbsRT <= cw.gibbsRT + 1.0e-8 && (!owGlobal ||
                (cc.stable && cw.stable && globalGap <= 1.0e-6));
            allPass = allPass && pass;
            gates << t << ',' << p << ',' << w << ',' << (pass ? "PASS" : "FAIL")
                << ',' << gap << ',' << globalGap
                << ",globally_stable_OW_requires_same_certified_branch_from_cold_and_continued_seeds\n";
            std::cout << "branch probe T=" << t << " P=" << p << " w=" << w
                << " cold-warm=" << gap << " status=" << (pass ? "PASS" : "FAIL") << std::endl;
        }
        std::cout << (allPass ? "BRANCH_CONSISTENCY_PASS" : "BRANCH_CONSISTENCY_BLOCKED") << std::endl;
        return allPass ? 0 : 2;
    } catch (const std::exception &e) {
        std::cerr << "branch probe failed: " << e.what() << std::endl;
        return 1;
    }
}
