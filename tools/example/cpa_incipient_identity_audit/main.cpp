// Diagnostic only. Reuse the frozen factory and production EOS/flash APIs.
#define main athabasca_reference_entrypoint
#include "../cpa_athabasca_bitumen_water/main.cpp"
#undef main

namespace {
const auto owMask()
{
    return MPMC::PhasePresence(static_cast<std::uint8_t>(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit));
}

double l1Gap(const Composition &a, const Composition &b)
{
    double value = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        value += std::abs(a[i] - b[i]);
    return value;
}
}

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3 && argc != 4)
            throw std::invalid_argument(
                "usage: identity_audit POINTWISE_CSV OUTPUT_DIR [FD_FACTOR]");
        const double fdFactor = argc == 4 ? std::stod(argv[3]) : 1.0e-5;
        if (!std::isfinite(fdFactor) || fdFactor <= 0.0 || fdFactor >= 0.01)
            throw std::invalid_argument("FD_FACTOR must be in (0,0.01)");
        std::ifstream in(argv[1]);
        if (!in) throw std::runtime_error("Cannot read pointwise CSV");
        std::string line;
        if (!std::getline(in, line)) throw std::runtime_error("Empty CSV");
        const auto header = splitCsv(line);
        const auto column = [&](const std::string &key) {
            const auto it = std::find(header.begin(), header.end(), key);
            if (it == header.end()) throw std::runtime_error("Missing " + key);
            return static_cast<std::size_t>(it - header.begin());
        };
        const auto tColumn = column("T_K");
        const auto pColumn = column("continuation_P_MPa");
        std::filesystem::create_directories(argv[2]);
        const std::filesystem::path out = argv[2];
        std::ofstream rows(out / "incipient_identity.csv");
        std::ofstream hess(out / "oil_tpd_hessian.csv");
        if (!rows || !hess) throw std::runtime_error("Cannot create outputs");
        rows << std::setprecision(17)
             << "T_K,boundary_P_MPa,offset_MPa,P_MPa,ow_converged,stability_valid,"
                "gas_unstable,trial_sum,candidate_xH2O,oil_xH2O,water_xH2O,"
                "candidate_oil_L1,candidate_water_L1,candidate_Z,vapor_Z,liquid_Z,"
                "oil_Z,water_Z,candidate_to_oil_density_ratio,"
                "candidate_to_water_density_ratio,candidate_fugacity_residual,"
                "candidate_tpd_RT\n";
        hess << std::setprecision(17)
             << "T_K,offset_MPa,row,column,H_forward,H_central\n";
        Eos eos = makeJiaCase1Cpa();
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = 0;
        options.maximumIterations = 240;
        options.maximumStabilityIterations = 160;
        options.cpaSelectGibbsMinimumRoot = true;
        Flash flash(eos, options);
        // Input must be the Amani-feed pointwise-probe CSV, not a fitted curve.
        const auto z = feedFromWaterMassFraction(0.441);
        int records = 0, failures = 0, temperatures = 0;
        while (std::getline(in, line))
        {
            if (line.empty()) continue;
            const auto fields = splitCsv(line);
            const double t = std::stod(fields.at(tColumn));
            const double pb = std::stod(fields.at(pColumn));
            if (!std::isfinite(t) || !std::isfinite(pb) || t <= 0 || pb <= 0 || pb >= 29.9)
                throw std::runtime_error("Invalid boundary P/T");
            ++temperatures;
            double walkP = 30.0;
            auto result = flash.flashRestricted(walkP * 1e6, t, z, owMask());
            const std::array<double, 5> offsets{{-0.020, 0.0002, 0.001, 0.005, 0.020}};
            for (double offset : offsets)
            {
                const double p = pb - offset;
                while (result.converged && walkP > p)
                {
                    walkP = std::max(p, walkP - 0.05);
                    result = flash.flashRestricted(
                        walkP * 1e6, t, z, owMask(), result.composition);
                }
                if (!result.converged || result.presence.bits() != owMask().bits())
                {
                    rows << t << ',' << pb << ',' << offset << ',' << p << ",0,0,0";
                    for (int k = 0; k < 15; ++k) rows << ",nan";
                    rows << '\n';
                    ++failures;
                    continue;
                }
                const auto st = flash.stabilityTest(
                    p * 1e6, t, z, result.presence, result.composition);
                const auto x = st.incipientComposition[1];
                double sum = 0.0;
                bool finite = true;
                for (double xi : x)
                {
                    sum += xi;
                    finite = finite && std::isfinite(xi) && xi > 0.0;
                }
                if (!st.valid || !finite || std::abs(sum - 1.0) >= 1e-8)
                {
                    ++failures;
                    continue;
                }
                const auto oil = eos.phaseResult(
                    p * 1e6, t, result.composition[0], MPMC::CompositionalPhase::Oil, false);
                const auto aq = eos.phaseResult(
                    p * 1e6, t, result.composition[2], MPMC::CompositionalPhase::Water, false);
                const auto selected = eos.phaseResult(
                    p * 1e6, t, x, MPMC::CompositionalPhase::Gas, true);
                const auto vapor = eos.phaseResult(
                    p * 1e6, t, x, MPMC::CompositionalPhase::Gas, false);
                const auto liquid = eos.phaseResult(
                    p * 1e6, t, x, MPMC::CompositionalPhase::Oil, false);
                double residual = 0.0, tpd = 0.0;
                for (std::size_t k = 0; k < 5; ++k)
                {
                    const double d = std::log(x[k] * selected.fugacityCoefficient[k])
                        - std::log(result.composition[0][k] * oil.fugacityCoefficient[k]);
                    residual = std::max(residual, std::abs(d + std::log(st.trialSum[1])));
                    tpd += x[k] * d;
                }
                rows << t << ',' << pb << ',' << offset << ',' << p << ",1,"
                     << st.valid << ',' << st.missingPhaseUnstable[1] << ','
                     << st.trialSum[1] << ',' << x[0] << ',' << result.composition[0][0]
                     << ',' << result.composition[2][0] << ',' << l1Gap(x, result.composition[0])
                     << ',' << l1Gap(x, result.composition[2]) << ',' << selected.compressibility
                     << ',' << vapor.compressibility << ',' << liquid.compressibility
                     << ',' << oil.compressibility << ',' << aq.compressibility
                     << ',' << oil.compressibility / selected.compressibility
                     << ',' << aq.compressibility / selected.compressibility
                     << ',' << residual << ',' << tpd << '\n';
                // Same-liquid-root finite differences of mu_i-mu_4 with
                // x_j += h, x_4 -= h. A negative symmetric-Hessian eigenvalue
                // diagnoses local composition instability, not vapor identity.
                const auto a = result.composition[0];
                const auto mu = [&](const Composition &c) {
                    const auto th = eos.phaseResult(
                        p * 1e6, t, c, MPMC::CompositionalPhase::Oil, false);
                    std::array<double, 4> values{};
                    for (int i = 0; i < 4; ++i)
                        values[i] = std::log(c[i] * th.fugacityCoefficient[i])
                            - std::log(c[4] * th.fugacityCoefficient[4]);
                    return values;
                };
                const auto base = mu(a);
                for (int j = 0; j < 4; ++j)
                {
                    const double h = fdFactor * std::min(a[j], a[4]);
                    auto plus = a, minus = a;
                    plus[j] += h; plus[4] -= h;
                    minus[j] -= h; minus[4] += h;
                    const auto mp = mu(plus), mm = mu(minus);
                    for (int i = 0; i < 4; ++i)
                        hess << t << ',' << offset << ',' << i << ',' << j << ','
                             << (mp[i] - base[i]) / h << ','
                             << (mp[i] - mm[i]) / (2.0 * h) << '\n';
                }
                ++records;
            }
            std::cout << "identity audit T=" << t << " completed" << std::endl;
        }
        std::cout << "records=" << records << " diagnostic_failures=" << failures << std::endl;
        // Completion is NOT a physical WLV-WL acceptance gate.
        return failures == 0 && temperatures > 0 && records == 5 * temperatures ? 0 : 2;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
