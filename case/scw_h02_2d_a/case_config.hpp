#pragma once
#include "../scw_kerogen_common/benchmark_common.hpp"
#include <indices/model_config.hpp>
#include <natural/thermo/thermodynamic_model.hpp>
#include <array>
#include <cstddef>

namespace H02A {
struct Model {
 static constexpr int numberOfComponents=1;
 static constexpr bool hasWater=true,hasWells=true,enableDissolution=false,
   enableAdsorption=false,enableLandTrapping=false;
 static constexpr auto phaseBehavior =
   MPMC::PhaseBehaviorModel::LegacyOilGasWithIndependentWater;
};
struct Grid {
 inline static constexpr const char *meshDirectory="";
 static constexpr int nx=60,ny=20,nz=1;
 static constexpr double lx=.30,ly=.10,lz=.010;
};
struct Rock {
 static constexpr double kx=1e-12,ky=1e-12,kz=1e-12,porosity=.25;
};
struct Fluid {
 static constexpr std::size_t N=1;
 static constexpr int heavyComponent=0;
 inline static constexpr std::array<const char*,1> componentNames{"OIL_HEAVY"};
 inline static constexpr std::array<double,1> criticalTemperature{982.864};
 inline static constexpr std::array<double,1> criticalPressure{.808714e6};
 inline static constexpr std::array<double,1> criticalVolume{1.793890e-3};
 inline static constexpr std::array<double,1> acentricFactor{1.215676};
 inline static constexpr std::array<double,1> molarMass{.660132};
 inline static constexpr std::array<std::array<double,1>,1> binaryInteraction{{{{0.0}}}};
 static constexpr auto thermodynamicModel=MPMC::CubicThermodynamicModel::PengRobinson;
 static constexpr int eosModelFlag=1;
 static constexpr double eosOmegaA=.45724,eosOmegaB=.07780;
 static constexpr double eosU=2.4142135623730951,eosW=-.4142135623730951;
 static constexpr double temperature=653.15;
 // [Oil, Gas, independent Water].  Water values are operating-condition
 // design closures; no component transfer between Water and Heavy exists.
 inline static constexpr std::array<double,3> surfaceDensity{768.2,768.2,364.5};
 inline static constexpr std::array<double,3> viscosity{.002,.002,5e-5};
 static constexpr double waterViscosity=5e-5;
 static constexpr double waterFormationVolumeFactor=1.0;
};
struct InitialState {
 static constexpr bool preserveReferenceValues=false;
 static constexpr double pressure=28e6,temperature=653.15;
 static constexpr double oilSaturation=1.0,gasSaturation=0.0,waterSaturation=0.0;
 inline static constexpr std::array<double,1> oilComposition{1.0};
 inline static constexpr std::array<double,1> gasComposition{1.0};
};
struct Dissolution {
 static constexpr int component=0;
 static constexpr double waterMolarMass=.01801528,salinityMolality=0,
   initialWaterCO2MoleFraction=0;
};
struct Land {static constexpr double constant=0;};
struct Adsorption {
 static constexpr double rockDensity=2650,standardPressure=101325,standardTemperature=288.15;
 inline static constexpr std::array<double,1> thetaMax{},coefficient{};
};
struct Numerics:ScwKerogen1D::Numerics {
 static constexpr double massResidualScale=2e-5,rateWellResidualFloor=2.5e-8;
};
struct Time:ScwKerogen1D::Time {
 static constexpr int numberOfSteps=200;
 static constexpr double dtDays=60.0/86400.0;
 static constexpr double minimumDtDays=1e-5/86400.0;
 static constexpr int maximumRetries=12;
 static constexpr double targetPVI=2.0;
};
struct Output:ScwKerogen1D::Output {
 static constexpr double producerEffectivePoreVolumeM3=7.5e-5;
 inline static constexpr const char *directory="./results";
 static constexpr std::size_t every=25;
};
struct Config {
 inline static constexpr const char *name="scw_h02_2d_a";
 using Model=H02A::Model;using Grid=H02A::Grid;using Rock=H02A::Rock;
 using Fluid=H02A::Fluid;using InitialState=H02A::InitialState;
 using Dissolution=H02A::Dissolution;using Land=H02A::Land;
 using Adsorption=H02A::Adsorption;using Numerics=H02A::Numerics;
 using Time=H02A::Time;using Output=H02A::Output;
};
}
