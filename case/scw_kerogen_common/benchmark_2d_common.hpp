#pragma once

#include "benchmark_common.hpp"

#include <array>

/**
 * @file benchmark_2d_common.hpp
 * @brief 60x20x1 旧数值 benchmark 的均质二维网格、岩石与井控。
 *
 * @warning 本文件的 1.20m x 0.10m x 0.10m、phi=0.35 与 1500/150 mD
 *          仅为历史 numerical benchmark。新的实验可复现 slab 只复用
 *          60x20x1 拓扑，物理尺寸/岩石参数必须来自实际装置和试件。
 *
 * The physical dimensions and pore volume are identical to the 60x1x1
 * control.  Refining only y therefore isolates transverse flow and the
 * point-well geometry from a change in injected pore volumes.
 */
namespace ScwKerogen2D
{

using ScwKerogen1D::bar;
using ScwKerogen1D::mD;
using ScwKerogen1D::secondsPerDay;
using ScwKerogen1D::temperature;
using ScwKerogen1D::initialPressure;
using ScwKerogen1D::waterCriticalTemperature;
using ScwKerogen1D::waterCriticalPressure;

struct Grid
{
    inline static constexpr const char *meshDirectory = "";
    static constexpr int nx = 60;
    static constexpr int ny = 20;
    static constexpr int nz = 1;
    static constexpr double lx = 1.20;
    static constexpr double ly = 0.10;
    static constexpr double lz = 0.10;
};

struct Rock
{
    static constexpr double kxValue = 1500.0 * mD;
    static constexpr double kyValue = 1500.0 * mD;
    static constexpr double kzValue = 150.0 * mD;
    static constexpr double porosityValue = 0.35;

    static constexpr double kx(int, int, int) noexcept { return kxValue; }
    static constexpr double ky(int, int, int) noexcept { return kyValue; }
    static constexpr double kz(int, int, int) noexcept { return kzValue; }
    static constexpr double porosity(int, int, int) noexcept
    {
        return porosityValue;
    }
};

inline constexpr double bulkVolume = Grid::lx * Grid::ly * Grid::lz;
inline constexpr double nominalPoreVolume = bulkVolume * Rock::porosityValue;
inline constexpr double poreVolumesPerDay = ScwKerogen1D::poreVolumesPerDay;
inline constexpr double reservoirRate =
    poreVolumesPerDay * nominalPoreVolume / secondsPerDay;

using Type = ScwKerogen1D::Type;
using Control = ScwKerogen1D::Control;
using InjectionPhase = ScwKerogen1D::InjectionPhase;
using Completion = ScwKerogen1D::Completion;
using WellDefinition = ScwKerogen1D::WellDefinition;

// With an even ny the geometric centre y=0.05 m lies between j=9 and j=10.
// The primary case uses the lower of the two mirror-equivalent cells.  Its
// centre is y=0.0475 m, only dy/2=2.5 mm from the exact short-edge centre.
inline constexpr int primaryCentreJ = Grid::ny / 2 - 1;
inline constexpr int mirrorCentreJ = Grid::ny / 2;
inline constexpr double centreOffset = 0.5 * Grid::ly / Grid::ny;

inline static constexpr std::array<WellDefinition, 2> wells{{
    {0, "SCW_INJ", Type::Injector, Control::ReservoirTotalRate,
     reservoirRate, initialPressure + 0.5 * bar,
     {0, primaryCentreJ, 0, 1}, 0.001, 0.0, InjectionPhase::Water, 0},
    {1, "PROD", Type::Producer, Control::ReservoirTotalRate,
     reservoirRate, initialPressure - 0.5 * bar,
     {-1, primaryCentreJ, 0, 1}, 0.001, 0.0, InjectionPhase::Oil, -1}
}};

using Numerics = ScwKerogen1D::Numerics;
using Time = ScwKerogen1D::Time;
using Output = ScwKerogen1D::Output;

} // namespace ScwKerogen2D
