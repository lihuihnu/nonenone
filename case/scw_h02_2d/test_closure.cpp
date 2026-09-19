#include "case_config.hpp"
#include "viscosity.hpp"
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <cassert>
#include <iostream>
using I=MPMC::ADIndices<MPMC::CompositionalModelConfig<2,true,true,false,false,false,MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>>;
int main(){auto f=MPMC::cases::makeFluidSystem<I,H02::Config>();using E=I::ValueType;f.flowViscosityOverride=[](E,const MPMC::FluidSystem<I>::Composition&w){return H02::viscosity(w[0],'C');};std::array<double,2>x{0.,1.};auto t=f.eos.phaseResult(28e6,653.15,x,MPMC::CompositionalPhase::Oil,false);auto [w,rho,mu]=f.eosFlowProperties(28e6,x,t.compressibility,653.15);assert(std::abs(mu-.002)<1e-14);assert(rho>0);assert(std::abs(H02::viscosity(1.,'C')-5e-5)<1e-14);assert(std::abs(H02::viscosity(.039,'B')-.002)<1e-14);std::cout<<"H02 production property hook PASS rhoH="<<rho<<" muH="<<mu<<std::endl;}
