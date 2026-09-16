/**
 * @file phase_role_canonicalization_test.cpp
 * @brief 单元测试：验证 `phase_role_canonicalization` 的核心语义、边界条件和回归行为。
 */
#include <common/units.hpp>
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/cubic_plus_association.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
namespace {
using Config=MPMC::CompositionalModelConfig<4,true,false,false,false,false,MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices=MPMC::ScalarIndices<Config>; using Eos=MPMC::CubicEquationOfState<Indices>; using Flash=MPMC::CubicThreePhaseFlash<Indices>; using C=std::array<double,4>;
constexpr std::array<double,4> tc{647.096,304.1282,190.564,723.0},pc{22.064e6,7.3773e6,4.5992e6,1.410e6},vc{5.6e-5,9.4e-5,9.9e-5,9.0e-4},om{0.3443,0.22394,0.01142,0.742},mw{0.01801528,0.0440095,0.016043,0.226441};
constexpr std::array<std::array<double,4>,4> kij{{{{0,.1896,.485,.5}},{{.1896,0,.12,.09}},{{.485,.12,0,0}},{{.5,.09,0,0}}}};
void req(bool q,const std::string&m){if(!q)throw std::runtime_error(m);} Eos pr(){MPMC::CompositionalMixture<Indices> m(tc,pc,vc,om,mw,kij); return Eos(.45724,.07780,std::move(m),5,2.414213562373095,-.414213562373095,1e-30);} Eos sw(){Eos e=pr(); Eos::SoreideWhitsonOptions o; o.waterComponent=0;o.salinityMolality=0;o.aqueousWaterBip[0]=[](double,double){return 0.;};o.aqueousWaterBip[1]=[](double T,double s){return MPMC::SoreideWhitsonCorrelations::co2AqueousBip(T,tc[1],s);};o.aqueousWaterBip[2]=[](double T,double s){return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(T,tc[2],om[2],s);};o.aqueousWaterBip[3]=[](double,double){return .5;};e.configureSoreideWhitson(std::move(o));return e;} Eos cpa(){Eos e=pr();Eos::CubicPlusAssociationOptions o;constexpr double R=MPMC::units::gasConstant;for(size_t i=0;i<4;++i){o.a0[i]=.42748*R*R*tc[i]*tc[i]/pc[i];o.b[i]=.08664*R*tc[i]/pc[i];o.c1[i]=.480+1.574*om[i]-.176*om[i]*om[i];}o.a0[0]=MPMC::StandardCpaWater4C::a0;o.b[0]=MPMC::StandardCpaWater4C::b;o.c1[0]=MPMC::StandardCpaWater4C::c1;o.associationEnergy[0]=MPMC::StandardCpaWater4C::epsilon;o.associationVolume[0]=MPMC::StandardCpaWater4C::beta;o.donorSites[0]=2;o.acceptorSites[0]=2;o.radialDistribution=MPMC::CpaRadialDistribution::Simplified;e.configureCubicPlusAssociation(std::move(o));return e;}
void check(const Eos&e,double T,double Pbar,const C&z,const char*name){MPMC::ThreePhaseFlashOptions o;o.waterComponent=0;Flash f(e,o);double P=Pbar*1e5;auto r=f.flash(P,T,z);req(r.converged&&r.presence.bits()==MPMC::PhasePresence::allBits,std::string(name)+": three-phase convergence");req(r.composition[2][0]>.99&&r.composition[2][0]>r.composition[0][0]&&r.composition[2][0]>r.composition[1][0],std::string(name)+": W identity");req(r.composition[1][1]>r.composition[0][1],std::string(name)+": gas must be richer in CO2 than oil");auto a=e.phaseResult(P,T,r.composition[0],MPMC::CompositionalPhase::Oil),b=e.phaseResult(P,T,r.composition[1],MPMC::CompositionalPhase::Gas),c=e.phaseResult(P,T,r.composition[2],MPMC::CompositionalPhase::Water);for(int i=0;i<4;++i){if(r.composition[0][i]<=1e-20||r.composition[1][i]<=1e-20||r.composition[2][i]<=1e-20)continue;double x=std::log(std::max(a.fugacity[i],1e-300)),y=std::log(std::max(b.fugacity[i],1e-300)),w=std::log(std::max(c.fugacity[i],1e-300));req(std::max({x,y,w})-std::min({x,y,w})<2e-7,std::string(name)+": fugacity closure");}}
}
int main(){try{check(sw(),342,49.370685698216974,{.55,.20,.20,.05},"SW G/W");check(cpa(),342,262.567915178188,{.10,.25,.50,.15},"CPA O/W");check(cpa(),280,49.370685698216974,{.05,.20,.20,.55},"CPA G/W");std::cout<<"Global O/G/W phase-role canonicalization: ALL PASS\n";return 0;}catch(const std::exception&e){std::cerr<<"Global O/G/W phase-role canonicalization: FAIL: "<<e.what()<<'\n';return 1;}}
