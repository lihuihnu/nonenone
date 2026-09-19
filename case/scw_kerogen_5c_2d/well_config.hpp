#pragma once
#include "case_config.hpp"
#include <limits>

namespace ScwKerogen5C2D {

struct WellDefinition : ScwKerogen1D::WellDefinition {
    double maximumBhp = std::numeric_limits<double>::quiet_NaN();
};

inline std::array<WellDefinition, 2> makeWellDefinitions(int nx, int ny)
{
    using namespace ScwKerogen1D;
    std::array<WellDefinition, 2> wells{};

    wells[0].id = 0;
    wells[0].name = "SCW_INJ";
    wells[0].type = Type::Injector;
    wells[0].control = Control::ReservoirTotalRate;
    wells[0].target = 2.5e-8;
    wells[0].initialBhp = 28.05e6;
    wells[0].maximumBhp = 30e6;
    wells[0].completion = {
        int(0.0125 / 0.30 * nx),
        int(0.0125 / 0.10 * ny),
        0, 1};
    wells[0].radius = 0.0002;
    wells[0].injectionPhase = InjectionPhase::Water;
    wells[0].injectedComponent = Fluid::waterComponent;

    wells[1].id = 1;
    wells[1].name = "PROD";
    wells[1].type = Type::Producer;
    wells[1].control = Control::Bhp;
    wells[1].target = 28e6;
    wells[1].initialBhp = 28e6;
    wells[1].completion = {
        int(0.2875 / 0.30 * nx),
        int(0.0875 / 0.10 * ny),
        0, 1};
    wells[1].radius = 0.0002;
    wells[1].injectedComponent = -1;

    return wells;
}

} // namespace ScwKerogen5C2D
