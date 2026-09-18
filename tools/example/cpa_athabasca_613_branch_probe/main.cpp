// Focused 613.2 K branch-map diagnostic for the Jia/Amani Figure-7 audit.
#define main athabasca_reference_entrypoint
#include "../cpa_athabasca_bitumen_water/main.cpp"
#undef main

namespace
{
const auto owPresence()
{
    return MPMC::PhasePresence(static_cast<std::uint8_t>(
        MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit));
}

double maxOwGap(const Flash::Result &a, const Flash::Result &b)
{
    if (!a.converged || !b.converged ||
        a.presence.bits() != owPresence().bits() ||
        b.presence.bits() != owPresence().bits())
        return std::numeric_limits<double>::infinity();
    double gap = 0.0;
    for (std::size_t phase : {std::size_t(0), std::size_t(2)})
        for (std::size_t i = 0; i < 5; ++i)
            gap = std::max(
                gap,
                std::abs(a.composition[phase][i] - b.composition[phase][i]));
    return gap;
}

double fixedRootGibbs(
    const Eos &eos,
    double pMPa,
    double t,
    const Flash::Result &r)
{
    if (!r.converged)
        return std::numeric_limits<double>::quiet_NaN();
    const double p = pMPa * 1.0e6;
    double g = 0.0;
    for (std::size_t k = 0; k < 3; ++k)
    {
        const auto role = static_cast<MPMC::CompositionalPhase>(k);
        if (!r.presence.contains(role) || !(r.phaseMoleFraction[k] > 0.0))
            continue;
        const auto thermo = eos.phaseResult(
            p, t, r.composition[k], role, false);
        for (std::size_t i = 0; i < 5; ++i)
        {
            const double x = r.composition[k][i];
            if (x > 0.0)
                g += r.phaseMoleFraction[k] * x *
                    (std::log(x) +
                     std::log(std::max(
                         thermo.fugacityCoefficient[i], 1.0e-300)));
        }
    }
    return g;
}

struct StabilityInfo
{
    bool valid{false};
    bool gasUnstable{false};
    double gasTrial{std::numeric_limits<double>::quiet_NaN()};
    double xOilWater{std::numeric_limits<double>::quiet_NaN()};
    double xWaterWater{std::numeric_limits<double>::quiet_NaN()};
};

StabilityInfo inspect(
    const Flash &flash,
    double pMPa,
    double t,
    const Composition &z,
    const Flash::Result &r)
{
    StabilityInfo info;
    if (!r.converged || r.presence.bits() != owPresence().bits())
        return info;
    const auto s = flash.stabilityTest(
        pMPa * 1.0e6, t, z, r.presence, r.composition);
    if (!s.valid)
        return info;
    const std::size_t gas = static_cast<std::size_t>(
        MPMC::phaseIndex(MPMC::CompositionalPhase::Gas));
    info.valid = true;
    info.gasUnstable = s.missingPhaseUnstable[gas];
    info.gasTrial = s.trialSum[gas];
    info.xOilWater = r.composition[0][water];
    info.xWaterWater = r.composition[2][water];
    return info;
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 2)
            throw std::invalid_argument("usage: branch613_probe OUTPUT_DIR");
        const std::filesystem::path out = argv[1];
        std::filesystem::create_directories(out);
        std::ofstream rows(out / "branch613_map.csv");
        if (!rows)
            throw std::runtime_error("Cannot open branch613_map.csv");
        rows << std::scientific << std::setprecision(12)
             << "P_MPa,cold_valid,cold_gas_unstable,cold_gas_trial_sum,"
                "warm_valid,warm_gas_unstable,warm_gas_trial_sum,"
                "cold_warm_composition_gap,cold_gibbs_RT,warm_gibbs_RT,"
                "cold_xH2O_oil,cold_xH2O_water,"
                "warm_xH2O_oil,warm_xH2O_water,"
                "global_converged,global_phase_code,global_stable,"
                "global_beta_gas,global_gibbs_RT\n";

        Eos eos = makeJiaCase1Cpa();
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = 0;
        options.maximumIterations = 240;
        options.maximumStabilityIterations = 160;
        options.cpaSelectGibbsMinimumRoot = true;
        Flash flash(eos, options);
        const Composition z = feedFromWaterMassFraction(0.441);
        constexpr double t = 613.2;

        // Establish the high-pressure liquid branch well above both candidate
        // crossings, then continue down in 0.01 MPa increments.
        double warmP = 18.60;
        auto warm = flash.flashRestricted(
            warmP * 1.0e6, t, z, owPresence());
        if (!warm.converged)
            throw std::runtime_error(
                "Cannot establish 613 K high-pressure O+W anchor.");

        bool anyBranchSplit = false;
        for (int n = 0; n <= 40; ++n)
        {
            const double targetP = 18.50 - 0.01 * n;
            while (warmP > targetP + 1.0e-12)
            {
                const double nextP = std::max(targetP, warmP - 0.01);
                auto next = flash.flashRestricted(
                    nextP * 1.0e6, t, z, owPresence(), warm.composition);
                if (!next.converged)
                    break;
                warm = next;
                warmP = nextP;
            }

            const auto cold = flash.flashRestricted(
                targetP * 1.0e6, t, z, owPresence());
            const auto coldInfo = inspect(flash, targetP, t, z, cold);
            const auto warmInfo = inspect(flash, targetP, t, z, warm);
            const double gap = maxOwGap(cold, warm);
            if (std::isfinite(gap) && gap > 1.0e-6)
                anyBranchSplit = true;

            const auto global = flash.flash(targetP * 1.0e6, t, z);
            bool globalStable = false;
            if (global.converged)
            {
                const auto s = flash.stabilityTest(
                    targetP * 1.0e6, t, z,
                    global.presence, global.composition);
                globalStable = s.valid && s.stable &&
                    maxMaterialClosure(z, global) <= 1.0e-8;
            }
            const std::size_t gas = static_cast<std::size_t>(
                MPMC::phaseIndex(MPMC::CompositionalPhase::Gas));

            rows << targetP << ','
                 << coldInfo.valid << ',' << coldInfo.gasUnstable << ','
                 << coldInfo.gasTrial << ','
                 << warmInfo.valid << ',' << warmInfo.gasUnstable << ','
                 << warmInfo.gasTrial << ',' << gap << ','
                 << fixedRootGibbs(eos, targetP, t, cold) << ','
                 << fixedRootGibbs(eos, targetP, t, warm) << ','
                 << coldInfo.xOilWater << ',' << coldInfo.xWaterWater << ','
                 << warmInfo.xOilWater << ',' << warmInfo.xWaterWater << ','
                 << global.converged << ','
                 << (global.converged ? global.presence.bits() : 0) << ','
                 << globalStable << ','
                 << (global.converged ? global.phaseMoleFraction[gas] : 0.0)
                 << ',' << fixedRootGibbs(eos, targetP, t, global) << '\n';
        }

        std::cout << "BRANCH613_MAP_COMPLETE branch_split="
                  << anyBranchSplit << std::endl;
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "branch613 probe failed: " << e.what() << std::endl;
        return 1;
    }
}
