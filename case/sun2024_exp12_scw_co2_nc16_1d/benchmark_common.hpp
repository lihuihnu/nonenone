/**
 * @file benchmark_common.hpp
 * @brief Sun et al. (2024) Exp. 12 派生超临界水-CO2-nC16 砂管算例公共参数。
 */
#pragma once

#include <common/math.hpp>
#include <natural/fluid_system.hpp>

#include <array>
#include <cstddef>
#include <limits>

namespace Sun2024Exp12
{

inline constexpr double pi = 3.14159265358979323846;
inline constexpr double mD = 9.869233e-16;
inline constexpr double day = 86400.0;
inline constexpr double minute = 60.0;

struct Grid
{
    // The paper's numerical model used 1 cm axial cells.  A 48x1x1 reduction
    // retains that axial resolution and the circular sandpack bulk volume.
    static constexpr int nx = 48;
    static constexpr int ny = 1;
    static constexpr int nz = 1;
    static constexpr double lx = 0.48;
    static constexpr double diameter = 0.039;
    static constexpr double crossSectionArea = pi * diameter * diameter / 4.0;
    static constexpr double ly = diameter;
    static constexpr double lz = crossSectionArea / ly;
    inline static constexpr const char *meshDirectory = "";
};

struct Rock
{
    static constexpr double permeability = 2000.0 * mD;
    static constexpr double porosityValue = 0.39;
    static constexpr double kx(int, int, int) { return permeability; }
    static constexpr double ky(int, int, int) { return permeability; }
    static constexpr double kz(int, int, int) { return permeability; }
    static constexpr double porosity(int, int, int) { return porosityValue; }
};

inline constexpr double bulkVolume =
    Grid::lx * Grid::crossSectionArea;
inline constexpr double poreVolume = bulkVolume * Rock::porosityValue;

// Exp. 12 reports 10 mL/min water and 2 mL/min CO2.  ReservoirTotalRate is
// used because the Natural well contract defines it as in-situ phase volume.
inline constexpr double waterInjectionRate = 10.0e-6 / minute;
inline constexpr double co2InjectionRate = 2.0e-6 / minute;
inline constexpr double totalInjectionRate = waterInjectionRate + co2InjectionRate;
inline constexpr double waterInjectionVolumeFraction =
    waterInjectionRate / totalInjectionRate;
inline constexpr double co2InjectionVolumeFraction =
    co2InjectionRate / totalInjectionRate;

inline constexpr double initialPressure = 24.0e6;
inline constexpr double temperature = 673.15;
inline constexpr double producerBhp = initialPressure;
// At 673.15 K the IF97 B23 boundary is 24.2356 MPa.  This rate case needs only
// a few kPa pressure gradient; the guard keeps accepted/Newton states inside
// the implemented, physically correct low-density Region 2 closure.
inline constexpr double injectorMaximumBhp = 24.20e6;
inline constexpr double targetSupercriticalFluidSaturation = 0.065;
inline constexpr double targetOilSaturation = 0.935;
inline constexpr double targetInjectedPoreVolumes = 4.0;
inline constexpr double finalTimeSeconds =
    targetInjectedPoreVolumes * poreVolume / totalInjectionRate;

struct Completion
{
    int i{0};
    int j{0};
    int kBegin{0};
    int kCount{1};
};

enum class Type { Injector, Producer };
enum class Control { Bhp, TotalRate, OilRate, GasRate, WaterRate };
enum class InjectionPhase { Oil, Gas, Water };

struct BaseWellDefinition
{
    int id{0};
    const char *name{""};
    Type type{Type::Producer};
    Control control{Control::Bhp};
    double target{0.0};
    double initialBhp{initialPressure};
    Completion completion{};
    double radius{0.001};
    double skin{0.0};
    InjectionPhase injectionPhase{InjectionPhase::Gas};
    double minimumBhp{std::numeric_limits<double>::infinity()};
    double maximumBhp{std::numeric_limits<double>::infinity()};
};

inline static constexpr std::array<BaseWellDefinition, 2> baseWells{{
    {0, "SCW_CO2_INJ", Type::Injector, Control::TotalRate,
     totalInjectionRate, initialPressure, {0, 0, 0, 1}, 0.001, 0.0,
     InjectionPhase::Gas, std::numeric_limits<double>::infinity(),
     injectorMaximumBhp},
    {1, "PROD", Type::Producer, Control::Bhp,
     producerBhp, producerBhp, {Grid::nx - 1, 0, 0, 1}, 0.001, 0.0,
     InjectionPhase::Oil, std::numeric_limits<double>::infinity(),
     std::numeric_limits<double>::infinity()}
}};

template <class Indices>
void applyRelativePermeability(MPMC::FluidSystem<Indices> &fluid)
{
    using Eval = typename Indices::ValueType;
    constexpr double swc = 0.0;
    constexpr double sorw = 0.15;
    constexpr double sgc = 0.0;
    constexpr double sorg = 0.10;

    fluid.waterRelativePermeability = [=](Eval sw) {
        const double value = MPMC::scalarValue(sw);
        if (value <= swc) return Eval(0.0);
        if (value >= 1.0 - sorw) return Eval(0.35);
        const Eval se = (sw - swc) / (1.0 - swc - sorw);
        return Eval(0.35) * se * se;
    };
    fluid.gasRelativePermeability = [=](Eval sg) {
        const double value = MPMC::scalarValue(sg);
        if (value <= sgc) return Eval(0.0);
        if (value >= 1.0 - sorg - swc) return Eval(0.85);
        const Eval se = (sg - sgc) / (1.0 - sorg - swc - sgc);
        return Eval(0.85) * se * se;
    };
    fluid.threePhaseOilRelativePermeability = [=](Eval, Eval so, Eval) {
        const double value = MPMC::scalarValue(so);
        if (value <= sorg) return Eval(0.0);
        if (value >= 1.0 - swc) return Eval(0.90);
        const Eval se = (so - sorg) / (1.0 - swc - sorg);
        return Eval(0.90) * se * se;
    };
}

} // namespace Sun2024Exp12
