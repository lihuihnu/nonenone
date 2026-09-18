// Focused near-boundary recovery audit at 593/603.6 K.
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

Flash::Result continuedOw(
    const Flash &flash,
    double t,
    const Composition &z,
    double targetP)
{
    double p = 30.0;
    auto r = flash.flashRestricted(p * 1.0e6, t, z, owPresence());
    if (!r.converged)
        return r;
    while (p > targetP)
    {
        const double next = std::max(targetP, p - 0.05);
        auto n = flash.flashRestricted(
            next * 1.0e6, t, z, owPresence(), r.composition);
        if (!n.converged)
            return n;
        r = n;
        p = next;
    }
    return r;
}

double fixedG(
    const Eos &eos,
    double pMPa,
    double t,
    const Flash::Result &r)
{
    if (!r.converged)
        return std::numeric_limits<double>::quiet_NaN();
    double g = 0.0;
    for (std::size_t k = 0; k < 3; ++k)
    {
        const auto role=static_cast<MPMC::CompositionalPhase>(k);
        if (!r.presence.contains(role) || r.phaseMoleFraction[k] <= 0.0)
            continue;
        const auto th=eos.phaseResult(
            pMPa*1.0e6,t,r.composition[k],role,false);
        for(std::size_t i=0;i<5;++i)
        {
            const double x=r.composition[k][i];
            if(x>0.0)
                g += r.phaseMoleFraction[k]*x*
                    (std::log(x)+std::log(std::max(
                        th.fugacityCoefficient[i],1.0e-300)));
        }
    }
    return g;
}
}

int main(int argc,char **argv)
{
    try
    {
        if(argc!=2)
            throw std::invalid_argument("usage: highT_onset_probe OUTPUT_DIR");
        const std::filesystem::path out=argv[1];
        std::filesystem::create_directories(out);
        std::ofstream rows(out/"highT_onset_recovery.csv");
        if(!rows) throw std::runtime_error("cannot create output");
        rows<<std::scientific<<std::setprecision(12)
            <<"T_K,boundary_P_MPa,offset_MPa,P_MPa,ow_converged,"
              "ow_stability_valid,ow_gas_unstable,ow_gas_trial_sum,"
              "global_converged,global_phase_code,global_stable,"
              "global_beta_gas,global_material_closure,"
              "incipient_xH2O,incipient_oil_L1,incipient_water_L1,"
              "incipient_selected_Z,incipient_vapor_Z,incipient_liquid_Z,"
              "global_gas_incipient_L1,global_gas_Z,global_gas_vapor_Z,"
              "global_gibbs_RT,ow_gibbs_RT,global_minus_ow_gibbs_RT\n";

        Eos eos=makeJiaCase1Cpa();
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent=0;
        options.maximumIterations=240;
        options.maximumStabilityIterations=160;
        options.cpaSelectGibbsMinimumRoot=true;
        Flash flash(eos,options);
        const Composition z=feedFromWaterMassFraction(0.441);
        const std::array<double,2> temps{{593.0,603.6}};
        const std::array<double,7> offsets{{
            0.005,0.010,0.020,0.050,0.100,0.200,0.400}};

        bool anyCertifiedThreeAtEachT=true;
        for(double t:temps)
        {
            const auto b=findWlvWlBoundaryContinuation(flash,t,z);
            if(!b.found) throw std::runtime_error("boundary not found");
            bool foundThree=false;
            for(double offset:offsets)
            {
                const double p=b.pressureMPa-offset;
                const auto ow=continuedOw(flash,t,z,p);
                bool owValid=false,gasUnstable=false;
                double trial=std::numeric_limits<double>::quiet_NaN();
                Composition incipient{};
                bool incipientValid=false;
                double incipientOilL1=std::numeric_limits<double>::quiet_NaN();
                double incipientWaterL1=std::numeric_limits<double>::quiet_NaN();
                double incipientSelectedZ=std::numeric_limits<double>::quiet_NaN();
                double incipientVaporZ=std::numeric_limits<double>::quiet_NaN();
                double incipientLiquidZ=std::numeric_limits<double>::quiet_NaN();
                if(ow.converged)
                {
                    const auto s=flash.stabilityTest(
                        p*1.0e6,t,z,ow.presence,ow.composition);
                    owValid=s.valid;
                    if(s.valid)
                    {
                        const std::size_t gas=static_cast<std::size_t>(
                            MPMC::phaseIndex(MPMC::CompositionalPhase::Gas));
                        gasUnstable=s.missingPhaseUnstable[gas];
                        trial=s.trialSum[gas];
                        incipient=s.incipientComposition[gas];
                        double sum=0.0;
                        incipientValid=true;
                        for(double x:incipient)
                        {
                            incipientValid=incipientValid&&std::isfinite(x)&&x>=0.0;
                            sum+=x;
                        }
                        incipientValid=incipientValid&&std::abs(sum-1.0)<=1.0e-8;
                        if(incipientValid)
                        {
                            incipientOilL1=0.0;
                            incipientWaterL1=0.0;
                            for(std::size_t i=0;i<5;++i)
                            {
                                incipientOilL1+=std::abs(
                                    incipient[i]-ow.composition[0][i]);
                                incipientWaterL1+=std::abs(
                                    incipient[i]-ow.composition[2][i]);
                            }
                            const auto selected=eos.phaseResult(
                                p*1.0e6,t,incipient,
                                MPMC::CompositionalPhase::Gas,true);
                            const auto vapor=eos.phaseResult(
                                p*1.0e6,t,incipient,
                                MPMC::CompositionalPhase::Gas,false);
                            const auto liquid=eos.phaseResult(
                                p*1.0e6,t,incipient,
                                MPMC::CompositionalPhase::Oil,false);
                            incipientSelectedZ=selected.compressibility;
                            incipientVaporZ=vapor.compressibility;
                            incipientLiquidZ=liquid.compressibility;
                        }
                    }
                }

                const auto global=flash.flash(p*1.0e6,t,z);
                bool globalStable=false;
                double closure=std::numeric_limits<double>::quiet_NaN();
                double betaGas=0.0;
                if(global.converged)
                {
                    closure=maxMaterialClosure(z,global);
                    const auto s=flash.stabilityTest(
                        p*1.0e6,t,z,global.presence,global.composition);
                    globalStable=s.valid&&s.stable&&closure<=1.0e-8;
                    betaGas=global.phaseMoleFraction[
                        static_cast<std::size_t>(
                            MPMC::phaseIndex(MPMC::CompositionalPhase::Gas))];
                    if(global.presence.bits()==MPMC::PhasePresence::allBits &&
                       globalStable &&
                       betaGas>options.phaseFractionTolerance)
                        foundThree=true;
                }
                double globalGasIncipientL1=
                    std::numeric_limits<double>::quiet_NaN();
                double globalGasZ=std::numeric_limits<double>::quiet_NaN();
                double globalGasVaporZ=std::numeric_limits<double>::quiet_NaN();
                if(global.converged &&
                   global.presence.contains(MPMC::CompositionalPhase::Gas))
                {
                    const std::size_t gas=static_cast<std::size_t>(
                        MPMC::phaseIndex(MPMC::CompositionalPhase::Gas));
                    globalGasZ=global.compressibility[gas];
                    const auto vapor=eos.phaseResult(
                        p*1.0e6,t,global.composition[gas],
                        MPMC::CompositionalPhase::Gas,false);
                    globalGasVaporZ=vapor.compressibility;
                    if(incipientValid)
                    {
                        globalGasIncipientL1=0.0;
                        for(std::size_t i=0;i<5;++i)
                            globalGasIncipientL1+=std::abs(
                                global.composition[gas][i]-incipient[i]);
                    }
                }
                const double gg=fixedG(eos,p,t,global);
                const double go=fixedG(eos,p,t,ow);
                rows<<t<<','<<b.pressureMPa<<','<<offset<<','<<p<<','
                    <<ow.converged<<','<<owValid<<','<<gasUnstable<<','<<trial<<','
                    <<global.converged<<','
                    <<(global.converged?global.presence.bits():0)<<','
                    <<globalStable<<','<<betaGas<<','<<closure<<','
                    <<(incipientValid?incipient[0]:
                       std::numeric_limits<double>::quiet_NaN())<<','
                    <<incipientOilL1<<','<<incipientWaterL1<<','
                    <<incipientSelectedZ<<','<<incipientVaporZ<<','
                    <<incipientLiquidZ<<','<<globalGasIncipientL1<<','
                    <<globalGasZ<<','<<globalGasVaporZ<<','
                    <<gg<<','<<go<<','<<(gg-go)<<'\n';
            }
            anyCertifiedThreeAtEachT =
                anyCertifiedThreeAtEachT && foundThree;
        }
        std::cout<<(anyCertifiedThreeAtEachT?
            "HIGH_T_ONSET_RECOVERY_CERTIFIED":
            "HIGH_T_ONSET_RECOVERY_BLOCKED")<<std::endl;
        return anyCertifiedThreeAtEachT?0:2;
    }
    catch(const std::exception &e)
    {
        std::cerr<<"highT onset probe failed: "<<e.what()<<std::endl;
        return 1;
    }
}
