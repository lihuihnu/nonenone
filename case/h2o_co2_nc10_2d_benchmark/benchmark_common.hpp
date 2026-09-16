/**
 * @file benchmark_common.hpp
 * @brief H2O-CO2-nC10 二维对比算例共享的网格、岩石与运行参数。
 */
#pragma once

#include <common/math.hpp>
#include <natural/fluid_system.hpp>

#include <array>
#include <cstddef>
#include <limits>

namespace BenchmarkCommon
{

inline constexpr double bar = 1.0e5;
inline constexpr double mD = 9.869233e-16;
inline constexpr double day = 86400.0;
inline constexpr double year = 365.25 * day;

struct Grid
{
    static constexpr int nx = 60;
    static constexpr int ny = 20;
    static constexpr int nz = 1;
    static constexpr double lx = 300.0;
    static constexpr double ly = 100.0;
    static constexpr double lz = 5.0;
    inline static constexpr const char *meshDirectory = "";
};

struct Rock
{
    static constexpr double permeability = 100.0 * mD;
    static constexpr double porosityValue = 0.20;
    static constexpr double kx(int, int, int) { return permeability; }
    static constexpr double ky(int, int, int) { return permeability; }
    static constexpr double kz(int, int, int) { return permeability; }
    static constexpr double porosity(int, int, int) { return porosityValue; }
};

inline constexpr double poreVolume = Grid::lx * Grid::ly * Grid::lz * Rock::porosityValue;
inline constexpr double injectionRate = 0.10 * poreVolume / year;
#ifndef MPMC_H2O_CO2_NC10_INITIAL_PRESSURE_PA
#define MPMC_H2O_CO2_NC10_INITIAL_PRESSURE_PA 5.16e6
#endif
#ifndef MPMC_H2O_CO2_NC10_TEMPERATURE_K
#define MPMC_H2O_CO2_NC10_TEMPERATURE_K 333.15
#endif
inline constexpr double initialPressure = MPMC_H2O_CO2_NC10_INITIAL_PRESSURE_PA;
inline constexpr double temperature = MPMC_H2O_CO2_NC10_TEMPERATURE_K;
inline constexpr double producerBhp = initialPressure - 0.50e6;
inline constexpr double injectorMaximumBhp = initialPressure + 2.0e6;
inline constexpr double targetWaterSaturation = 0.20;
inline constexpr double targetOilSaturation = 0.80;

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
    double radius{0.10};
    double skin{0.0};
    InjectionPhase injectionPhase{InjectionPhase::Gas};
    double minimumBhp{std::numeric_limits<double>::infinity()};
    double maximumBhp{std::numeric_limits<double>::infinity()};
};

inline static constexpr std::array<BaseWellDefinition, 2> baseWells{{
    {0, "CO2_INJ", Type::Injector, Control::TotalRate,
     injectionRate, initialPressure, {0, 9, 0, 1}, 0.10, 0.0,
     InjectionPhase::Gas, std::numeric_limits<double>::infinity(), injectorMaximumBhp},
    {1, "PROD", Type::Producer, Control::Bhp,
     producerBhp, producerBhp, {59, 9, 0, 1}, 0.10, 0.0,
     InjectionPhase::Oil, std::numeric_limits<double>::infinity(),
     std::numeric_limits<double>::infinity()}
}};

template <class Indices>
void applyCommonRelativePermeability(MPMC::FluidSystem<Indices> &fluid)
{
    using Eval = typename Indices::ValueType;
    constexpr double swc = 0.10;
    constexpr double sorw = 0.15;
    constexpr double sgc = 0.02;
    constexpr double sorg = 0.10;

    fluid.waterRelativePermeability = [=](Eval sw) {
        const double value = MPMC::scalarValue(sw);
        if (value <= swc) return Eval(0.0);
        if (value >= 1.0 - sorw) return Eval(0.30);
        const Eval se = (sw - swc) / (1.0 - swc - sorw);
        return Eval(0.30) * se * se;
    };
    fluid.gasRelativePermeability = [=](Eval sg) {
        const double value = MPMC::scalarValue(sg);
        if (value <= sgc) return Eval(0.0);
        if (value >= 1.0 - sorg - swc) return Eval(0.80);
        const Eval se = (sg - sgc) / (1.0 - sorg - swc - sgc);
        return Eval(0.80) * se * se;
    };
    fluid.threePhaseOilRelativePermeability = [=](Eval, Eval so, Eval) {
        const double value = MPMC::scalarValue(so);
        if (value <= sorg) return Eval(0.0);
        if (value >= 1.0 - swc) return Eval(0.80);
        const Eval se = (so - sorg) / (1.0 - swc - sorg);
        return Eval(0.80) * se * se;
    };
}

} // namespace BenchmarkCommon
