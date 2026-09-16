/**
 * @file benchmark_common.hpp
 * @brief 三维超临界水-CO2-nC4-nC16 运移算例的网格、岩石与注采基线。
 */
#pragma once

#include <common/math.hpp>
#include <natural/fluid_system.hpp>

#include <array>
#include <cstddef>
#include <limits>

namespace ScwMigration3D
{

inline constexpr double mD = 9.869233e-16;
inline constexpr double day = 86400.0;
inline constexpr double year = 365.25 * day;
inline constexpr double waterCriticalTemperature = 647.096;
inline constexpr double waterCriticalPressure = 22.064e6;

struct Grid
{
    // 24 x 16 x 6 = 2304 cells; dx/dy/dz = 10/10/6 m.
    static constexpr int nx = 24;
    static constexpr int ny = 16;
    static constexpr int nz = 6;
    static constexpr double lx = 240.0;
    static constexpr double ly = 160.0;
    static constexpr double lz = 36.0;
    inline static constexpr const char *meshDirectory = "";
};

/**
 * @brief 分层砂岩：斜向高渗通道穿过上下层，中部隔层仅保留局部窗口。
 *
 * k=2 是低渗隔层；i=10..14 且靠近通道中心的位置形成导流窗口。该结构使
 * 下部注入流体先沿通道推进，再由窗口和重力共同控制向上部生产层的迁移。
 */
struct Rock
{
    inline static constexpr std::array<double, Grid::nz> layerKxMd{
        700.0, 550.0, 12.0, 400.0, 280.0, 180.0};
    inline static constexpr std::array<double, Grid::nz> layerKyMd{
        520.0, 410.0, 9.0, 300.0, 210.0, 135.0};
    inline static constexpr std::array<double, Grid::nz> layerKzMd{
        55.0, 42.0, 0.30, 30.0, 20.0, 12.0};
    inline static constexpr std::array<double, Grid::nz> layerPorosity{
        0.235, 0.225, 0.120, 0.215, 0.200, 0.185};

    static constexpr int channelCenter(int i)
    {
        return 2 + (11 * i) / (Grid::nx - 1);
    }

    static constexpr bool inChannel(int i, int j, int k)
    {
        const int center = channelCenter(i);
        const int offset = j > center ? j - center : center - j;
        return k != 2 && offset <= 1;
    }

    static constexpr bool inBaffleWindow(int i, int j, int k)
    {
        if (k != 2 || i < 10 || i > 14)
            return false;
        const int center = channelCenter(i);
        const int offset = j > center ? j - center : center - j;
        return offset <= 1;
    }

    static constexpr double kx(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 180.0 * mD;
        const double multiplier = inChannel(i, j, k) ? 2.0 : 1.0;
        return layerKxMd[static_cast<std::size_t>(k)] * multiplier * mD;
    }

    static constexpr double ky(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 135.0 * mD;
        const double multiplier = inChannel(i, j, k) ? 2.0 : 1.0;
        return layerKyMd[static_cast<std::size_t>(k)] * multiplier * mD;
    }

    static constexpr double kz(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 14.0 * mD;
        const double multiplier = inChannel(i, j, k) ? 1.5 : 1.0;
        return layerKzMd[static_cast<std::size_t>(k)] * multiplier * mD;
    }

    static constexpr double porosity(int i, int j, int k)
    {
        if (inBaffleWindow(i, j, k))
            return 0.19;
        const double channelIncrement = inChannel(i, j, k) ? 0.012 : 0.0;
        return layerPorosity[static_cast<std::size_t>(k)] + channelIncrement;
    }
};

inline constexpr double nominalPorosity = 0.20;
inline constexpr double nominalPoreVolume =
    Grid::lx * Grid::ly * Grid::lz * nominalPorosity;
inline constexpr double injectionPoreVolumesPerYear = 0.025;
inline constexpr double totalInjectionRate =
    injectionPoreVolumesPerYear * nominalPoreVolume / year;
inline constexpr double targetInjectedPoreVolumes = 0.25;

inline constexpr double temperature = 673.15;
inline constexpr double initialPressure = 25.0e6;
// Neutral initial guesses and safety limits for the equal-rate well pair.
inline constexpr double producerInitialBhp = 24.5e6;
inline constexpr double producerMinimumBhp = 15.0e6;
inline constexpr double injectorMaximumBhp = 35.0e6;

inline constexpr double initialWaterSaturation = 0.15;
inline constexpr double initialOilSaturation = 0.85;
inline constexpr double initialGasSaturation = 0.0;
inline constexpr double defaultScwInjectionVolumeFraction = 0.75;

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
    {0, "SCW_CO2_INJ", Type::Injector, Control::TotalRate,
     totalInjectionRate, initialPressure, {2, Rock::channelCenter(2), 0, 2},
     0.10, 0.0, InjectionPhase::Gas,
     std::numeric_limits<double>::infinity(), injectorMaximumBhp},
    {1, "PROD", Type::Producer, Control::TotalRate,
     totalInjectionRate, producerInitialBhp,
     {Grid::nx - 3, Rock::channelCenter(Grid::nx - 3), 3, 3},
     0.10, 0.0, InjectionPhase::Oil,
     producerMinimumBhp,
     std::numeric_limits<double>::infinity()}
}};

template <class Indices>
void applyRelativePermeability(MPMC::FluidSystem<Indices> &fluid)
{
    using Eval = typename Indices::ValueType;
    constexpr double swc = 0.10;
    constexpr double sorw = 0.15;
    constexpr double sgc = 0.01;
    constexpr double sorg = 0.12;

    fluid.waterRelativePermeability = [=](Eval sw) {
        const double value = MPMC::scalarValue(sw);
        if (value <= swc) return Eval(0.0);
        if (value >= 1.0 - sorw) return Eval(0.45);
        const Eval se = (sw - swc) / (1.0 - swc - sorw);
        return Eval(0.45) * se * se;
    };
    fluid.gasRelativePermeability = [=](Eval sg) {
        const double value = MPMC::scalarValue(sg);
        if (value <= sgc) return Eval(0.0);
        if (value >= 1.0 - sorg - swc) return Eval(0.90);
        const Eval se = (sg - sgc) / (1.0 - sorg - swc - sgc);
        return Eval(0.90) * se * se;
    };
    fluid.threePhaseOilRelativePermeability = [=](Eval, Eval so, Eval) {
        const double value = MPMC::scalarValue(so);
        if (value <= sorg) return Eval(0.0);
        if (value >= 1.0 - swc) return Eval(0.85);
        const Eval se = (so - sorg) / (1.0 - swc - sorg);
        return Eval(0.85) * se * se;
    };
}

} // namespace ScwMigration3D
