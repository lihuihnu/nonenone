#pragma once
#include "../scw_kerogen_common/benchmark_common.hpp"
#include <indices/model_config.hpp>
#include <natural/thermo/thermodynamic_model.hpp>
#include <array>
namespace H02 {
struct Model {
 static constexpr int numberOfComponents=2;
 static constexpr bool hasWater=true,hasWells=true,enableDissolution=false,enableAdsorption=false,enableLandTrapping=false;
 static constexpr auto phaseBehavior=MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase;
};
struct Grid {
 inline static constexpr const char *meshDirectory="";
 static constexpr int nx=60,ny=20,nz=1;
 static constexpr double lx=.30,ly=.10,lz=.010;
};
struct Rock {static constexpr double kx=1e-12,ky=1e-12,kz=1e-12,porosity=.25;};
struct Fluid {
 static constexpr int waterComponent=0,heavyComponent=1;
 static constexpr std::size_t N=2;
 inline static constexpr std::array<const char*,2> componentNames{"H2O","OIL_HEAVY"};
 inline static constexpr std::array<double,2> criticalTemperature{647.096,982.864};
 inline static constexpr std::array<double,2> criticalPressure{22.064e6,.808714e6};
 inline static constexpr std::array<double,2> criticalVolume{7.49587808839913e-5,1.793890e-3};
 inline static constexpr std::array<double,2> acentricFactor{.344,1.215676};
 inline static constexpr std::array<double,2> molarMass{.01801528,.660132};
 inline static constexpr std::array<std::array<double,2>,2> binaryInteraction{{{{0,.30}},{{.30,0}}}};
 static constexpr auto thermodynamicModel=MPMC::CubicThermodynamicModel::PengRobinson;
 static constexpr int eosModelFlag=1;
 static constexpr double eosOmegaA=.45724,eosOmegaB=.07780,eosU=2.4142135623730951,eosW=-.4142135623730951;
 static constexpr double temperature=653.15;
 inline static constexpr std::array<double,3> surfaceDensity{768.2,364.5,364.5},viscosity{.002,5e-5,5e-5};
 static constexpr double waterViscosity=5e-5,waterFormationVolumeFactor=1.;
};
struct InitialState {static constexpr double pressure=28e6,temperature=653.15; inline static constexpr std::array<double,2> overallComposition{0,1};};
struct Dissolution {static constexpr int component=1;static constexpr double waterMolarMass=.01801528,salinityMolality=0,initialWaterCO2MoleFraction=0;};
struct Land {static constexpr double constant=0;};
struct Adsorption {static constexpr double rockDensity=2650,standardPressure=101325,standardTemperature=288.15;inline static constexpr std::array<double,2> thetaMax{},coefficient{};};
struct Numerics:ScwKerogen1D::Numerics {static constexpr double massResidualScale=2e-5,rateWellResidualFloor=2.5e-8;};
struct Time:ScwKerogen1D::Time {static constexpr int numberOfSteps=2;static constexpr double dtDays=.1/86400.,minimumDtDays=1e-5/86400.;static constexpr int maximumRetries=12;};
struct Output:ScwKerogen1D::Output {static constexpr double producerEffectivePoreVolumeM3=7.5e-5;inline static constexpr const char *directory="./results";};
struct Config {inline static constexpr const char *name="scw_h02_2d";using Model=H02::Model;using Grid=H02::Grid;using Rock=H02::Rock;using Fluid=H02::Fluid;using InitialState=H02::InitialState;using Dissolution=H02::Dissolution;using Land=H02::Land;using Adsorption=H02::Adsorption;using Numerics=H02::Numerics;using Time=H02::Time;using Output=H02::Output;};
}
