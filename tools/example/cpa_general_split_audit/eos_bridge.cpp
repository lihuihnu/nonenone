// Diagnostic C ABI over the unchanged production CPA factory. No phase-slot cap.
#define main athabasca_benchmark_main
#include "../cpa_athabasca_bitumen_water/main.cpp"
#undef main
#include <cstring>
namespace {
Eos& frozenEos() { static Eos eos = makeJiaCase1Cpa(); return eos; }
thread_local std::string lastError;
}
extern "C" {
const char* split_last_error() { return lastError.c_str(); }
int split_mu(double t, double p_mpa, int vapor, const double* x,
             double* mu, double* zfactor) noexcept {
    try {
        Composition c{}; double sum=0;
        for (int i=0;i<5;++i) { if (!(x[i]>0)||!std::isfinite(x[i])) throw std::invalid_argument("Positive composition required"); c[i]=x[i]; sum+=x[i]; }
        if(std::abs(sum-1)>1e-8) throw std::invalid_argument("Normalized composition required");
        const auto r=frozenEos().phaseResult(p_mpa*1e6,t,c,
            vapor?MPMC::CompositionalPhase::Gas:MPMC::CompositionalPhase::Oil,false);
        for(int i=0;i<5;++i) {mu[i]=std::log(c[i])+std::log(r.fugacityCoefficient[i]); if(!std::isfinite(mu[i]))throw std::runtime_error("Nonfinite chemical potential");}
        *zfactor=r.compressibility; return 0;
    } catch(const std::exception&e) {lastError=e.what();return 1;}
}
int split_feed(double* z) noexcept { const auto c=feedFromWaterMassFraction(0.441);std::copy(c.begin(),c.end(),z);return 0; }
int split_original_flash(double t,double p_mpa,double* n,int* roots,int* missing_stable) noexcept {
    try {
        MPMC::ThreePhaseFlashOptions o;o.waterComponent=0;o.maximumIterations=240;o.maximumStabilityIterations=160;o.cpaSelectGibbsMinimumRoot=true;
        Flash f(frozenEos(),o);auto z=feedFromWaterMassFraction(0.441);auto r=f.flash(p_mpa*1e6,t,z);
        if(!r.converged)throw std::runtime_error("Original three-slot flash did not converge");
        const auto s=f.stabilityTest(p_mpa*1e6,t,z,r.presence,r.composition);*missing_stable=s.valid&&s.stable;
        int q=0;for(int k=0;k<3;++k)if(r.presence.contains(static_cast<MPMC::CompositionalPhase>(k))){roots[q]=(k==1);for(int i=0;i<5;++i)n[5*q+i]=r.phaseMoleFraction[k]*r.composition[k][i];++q;}return q;
    } catch(const std::exception&e){lastError=e.what();return -1;}
}
}
