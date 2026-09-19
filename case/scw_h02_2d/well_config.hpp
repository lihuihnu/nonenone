#pragma once
#include "case_config.hpp"
#include <limits>
namespace H02 {
struct WellDefinition:ScwKerogen1D::WellDefinition {double maximumBhp=std::numeric_limits<double>::quiet_NaN();};
inline std::array<WellDefinition,2> makeWellDefinitions(int nx,int ny) {
 using namespace ScwKerogen1D;
 std::array<WellDefinition,2> w{};
 w[0].id=0;w[0].name="SCW_INJ";w[0].type=Type::Injector;w[0].control=Control::ReservoirTotalRate;w[0].target=2.5e-8;w[0].initialBhp=28.05e6;w[0].maximumBhp=30e6;w[0].completion={int(.0125/.30*nx),int(.0125/.10*ny),0,1};w[0].radius=.0002;w[0].injectionPhase=InjectionPhase::Water;w[0].injectedComponent=0;
 w[1].id=1;w[1].name="PROD";w[1].type=Type::Producer;w[1].control=Control::Bhp;w[1].target=28e6;w[1].initialBhp=28e6;w[1].completion={int(.2875/.30*nx),int(.0875/.10*ny),0,1};w[1].radius=.0002;w[1].injectedComponent=-1;
 return w;
}
}
