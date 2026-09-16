#pragma once

#include <array>

/**
 * @file benchmark_common.hpp
 * @brief 超临界水驱替干酪根裂解产物一维算例的共用网格与控制参数。
 */
namespace ScwKerogen1D
{

inline constexpr double bar = 1.0e5;
inline constexpr double mD = 9.869232667160130e-16;
inline constexpr double secondsPerDay = 86400.0;
inline constexpr double temperature = 653.2;
inline constexpr double initialPressure = 277.4 * bar;
inline constexpr double waterCriticalTemperature = 647.096;
inline constexpr double waterCriticalPressure = 22.064e6;

struct Grid
{
    inline static constexpr const char *meshDirectory = "";
    static constexpr int nx = 60;
    static constexpr int ny = 1;
    static constexpr int nz = 1;
    static constexpr double lx = 1.20;
    static constexpr double ly = 0.10;
    static constexpr double lz = 0.10;
};

struct Rock
{
    static constexpr double kx = 1500.0 * mD;
    static constexpr double ky = 1500.0 * mD;
    static constexpr double kz = 150.0 * mD;
    static constexpr double porosity = 0.35;
};

inline constexpr double bulkVolume = Grid::lx * Grid::ly * Grid::lz;
inline constexpr double nominalPoreVolume = bulkVolume * Rock::porosity;
inline constexpr double poreVolumesPerDay = 0.25;
inline constexpr double reservoirRate =
    poreVolumesPerDay * nominalPoreVolume / secondsPerDay;

enum class Type { Injector, Producer };
enum class Control { Bhp, TotalRate, OilRate, GasRate, WaterRate,
                     ReservoirTotalRate };
enum class InjectionPhase { Oil, Gas, Water };

struct Completion
{
    int i{0};
    int j{0};
    int kBegin{0};
    int kCount{-1};
};

struct WellDefinition
{
    int id{0};
    const char *name{""};
    Type type{Type::Producer};
    Control control{Control::Bhp};
    double target{0.0};
    double initialBhp{initialPressure};
    Completion completion{};
    double radius{0.004};
    double skin{0.0};
    InjectionPhase injectionPhase{InjectionPhase::Water};
    int injectedComponent{-1};
};

// H2O is component zero in both staged fluids.  Equal in-situ rates isolate
// compositional displacement from bulk pressure depletion.
inline static constexpr std::array<WellDefinition, 2> wells{{
    {0, "SCW_INJ", Type::Injector, Control::ReservoirTotalRate,
     reservoirRate, initialPressure + 0.5 * bar,
     {0, 0, 0, 1}, 0.004, 0.0, InjectionPhase::Water, 0},
    {1, "PROD", Type::Producer, Control::ReservoirTotalRate,
     reservoirRate, initialPressure - 0.5 * bar,
     {-1, 0, 0, 1}, 0.004, 0.0, InjectionPhase::Oil, -1}
}};

struct Numerics
{
    static constexpr double fugacityScalingFactor = 1.0;
    static constexpr bool useVariableBounds = false;
    // Characteristic scales are matched to this laboratory core and well
    // rate.  Without them the O(1e-6 kg/s) component equations and O(1e-8
    // m3/s) rate equations are numerically invisible beside closure rows.
    static constexpr double pressureScale = 1.0e7;
    static constexpr double compositionScale = 1.0;
    static constexpr double saturationScale = 1.0;
    static constexpr double massResidualScale = 1.0e-5;
    static constexpr double fugacityResidualScale = 1.0;
    static constexpr double closureResidualScale = 1.0;
    static constexpr double rateWellResidualFloor = reservoirRate;
    static constexpr bool enableSnesStagnationGuard = true;
    static constexpr int snesStagnationMinimumIterations = 8;
    static constexpr int snesStagnationWindow = 5;
    static constexpr double snesStagnationRelativeImprovement = 1.0e-4;
};

struct Time
{
    // 40 common output targets span 1 injected pore volume.
    static constexpr int numberOfSteps = 40;
    static constexpr double dtDays = 0.10;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = 1.0e-7;
    static constexpr double cutFactor = 0.5;
    static constexpr double growthFactor = 1.25;
    static constexpr double difficultShrinkFactor = 0.8;
    static constexpr int easyNonlinearIterations = 6;
    static constexpr int difficultNonlinearIterations = 14;
    static constexpr int maximumRetries = 16;
    static constexpr int maximumWellControlIterations = 8;
    static constexpr bool printAdaptiveSteps = true;
};

struct Output
{
    inline static constexpr const char *directory = "./results";
    static constexpr std::size_t every = 1;
    static constexpr bool printWells = true;
    static constexpr bool printWellPhaseDetails = true;
    static constexpr bool writeWellHistory = true;
    static constexpr bool printInventory = true;
    static constexpr bool enableComponentMassBalance = true;
    static constexpr bool printComponentMassBalance = true;
    static constexpr bool writeComponentMassBalance = true;
    static constexpr bool writeInventoryHistory = true;
    static constexpr bool writeMassTotals = true;
    static constexpr bool writeSolutionSnapshots = true;
    static constexpr bool writePhaseStateSnapshots = true;
    static constexpr bool printNewtonIterations = true;
    static constexpr bool writeSolverHistory = true;
    static constexpr bool printWellControlSwitches = true;
    static constexpr bool writeWellControlSwitchHistory = true;
    static constexpr bool printReservoirDiagnostics = true;
    static constexpr bool writeReservoirDiagnostics = true;
    static constexpr bool printFinalSummary = true;
    static constexpr bool writeFinalSummary = true;
};

} // namespace ScwKerogen1D
