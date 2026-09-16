/**
 * @file control_switch_output.hpp
 * @brief 井控制方式自动切换的独立日志输出。
 */
#pragma once

#include <common/console.hpp>

#include <well/control.hpp>
#include <well/state.hpp>
#include <well/types.hpp>

#include <petscsys.h>

#include <array>
#include <cstddef>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>

namespace MPMC
{

struct WellControlSwitchOutputOptions final
{
    std::filesystem::path resultDirectory{"./results"};
    std::filesystem::path historyFile{"well_control_switches.csv"};
    bool print{true};
    bool writeHistory{true};
    double secondsPerDisplayUnit{86400.0};
    std::string displayTimeUnit{"day"};
    double pressureScale{1.0e5};
    std::string pressureUnit{"bar"};
};

/**
 * @brief 自动井控制切换的专用输出。
 *
 * A switch is rare and important enough to deserve its own log line.  The
 * record contains the well identity, previous/new control, trigger, measured
 * value, limit, BHP and current phase rates.
 */
template <class Indices>
class WellControlSwitchOutput final
{
public:
    explicit WellControlSwitchOutput(
        WellControlSwitchOutputOptions options,
        PetscMPIInt rank)
        : options_(std::move(options)), rank_(rank)
    {
        if (rank_ != 0)
            return;

        std::filesystem::create_directories(options_.resultDirectory);
        if (!options_.writeHistory)
            return;

        const auto file = options_.resultDirectory / options_.historyFile;
        std::ofstream stream(file, std::ios::out | std::ios::trunc);
        if (!stream)
            throw std::runtime_error("Failed to create well-control switch history: " + file.string());

        stream
            << "time,well_id,well_name,well_type,old_control,old_target,new_control,new_target,"
            << "reason,measured,limit,bhp,q_oil_surface,q_gas_surface,q_water_surface\n";
    }

    [[nodiscard]] std::size_t switchCount() const noexcept
    {
        return switchCount_;
    }

    template <class Well>
    void write(
        double timeSeconds,
        const Well &well,
        const WellState<Indices> &state,
        const WellControlUpdate &update) const
    {
        if (rank_ != 0 || !update.changed)
            return;

        ++switchCount_;

        const bool pressureTrigger =
            update.reason == WellControlSwitchReason::MaximumBhp ||
            update.reason == WellControlSwitchReason::MinimumBhp;

        const double oldTarget = scaledTarget_(update.previousControl, update.previousTarget);
        const double newTarget = scaledTarget_(update.newControl, update.newTarget);
        const double measured = pressureTrigger
            ? update.measuredValue / options_.pressureScale
            : update.measuredValue;
        const double limit = pressureTrigger
            ? update.limitValue / options_.pressureScale
            : update.limitValue;

        if (options_.print)
        {
            const std::string triggerUnit =
                pressureTrigger ? options_.pressureUnit : "m3/s";
            const std::string transition =
                std::string(wellControlName(update.previousControl)) + "(" +
                consoleNumber(oldTarget) + ") -> " +
                std::string(wellControlName(update.newControl)) + "(" +
                consoleNumber(newTarget) + ")";
            const std::string rates =
                consoleScientific(phase_(state.surfacePhaseRate, Indices::Phase::liquid)) + " / " +
                consoleScientific(phase_(state.surfacePhaseRate, Indices::Phase::vapor)) + " / " +
                consoleScientific(water_(state.surfacePhaseRate));

            ConsoleSection section("WELL CONTROL SWITCH");
            section.row("Time", consoleNumber(timeSeconds / options_.secondsPerDisplayUnit),
                        options_.displayTimeUnit)
                .row("Well", std::to_string(well.id) + " / " + well.name + " / " +
                     std::string(wellTypeName(well.type)))
                .separator()
                .row("Control transition", transition)
                .row("Trigger", std::string(wellControlSwitchReasonName(update.reason)))
                .row("Measured / limit",
                     consoleNumber(measured) + " / " + consoleNumber(limit), triggerUnit)
                .row("BHP", consoleNumber(state.bottomHolePressure / options_.pressureScale),
                     options_.pressureUnit)
                .row("Surface q [O / G / W]", rates, "m3/s");
            const std::string text = section.str();
            PetscPrintf(PETSC_COMM_SELF, "%s", text.c_str());
        }

        if (options_.writeHistory)
            append_(timeSeconds, well, state, update, oldTarget, newTarget, measured, limit);
    }

private:
    [[nodiscard]] double scaledTarget_(WellControl control, double target) const noexcept
    {
        return control == WellControl::Bhp ? target / options_.pressureScale : target;
    }

    template <std::size_t N>
    [[nodiscard]] static double phase_(const std::array<double, N> &values, int phase)
    {
        return values[static_cast<std::size_t>(phase)];
    }

    template <std::size_t N>
    [[nodiscard]] static double water_(const std::array<double, N> &values)
    {
        if constexpr (Indices::hasWater)
            return phase_(values, Indices::Phase::water);
        return 0.0;
    }

    template <class Well>
    void append_(
        double timeSeconds,
        const Well &well,
        const WellState<Indices> &state,
        const WellControlUpdate &update,
        double oldTarget,
        double newTarget,
        double measured,
        double limit) const
    {
        const auto file = options_.resultDirectory / options_.historyFile;
        std::ofstream stream(file, std::ios::out | std::ios::app);
        if (!stream)
            throw std::runtime_error("Failed to append well-control switch history: " + file.string());

        stream << std::setprecision(16)
               << timeSeconds / options_.secondsPerDisplayUnit << ','
               << well.id << ',' << well.name << ',' << wellTypeName(well.type) << ','
               << wellControlName(update.previousControl) << ',' << oldTarget << ','
               << wellControlName(update.newControl) << ',' << newTarget << ','
               << wellControlSwitchReasonName(update.reason) << ','
               << measured << ',' << limit << ','
               << state.bottomHolePressure / options_.pressureScale << ','
               << phase_(state.surfacePhaseRate, Indices::Phase::liquid) << ','
               << phase_(state.surfacePhaseRate, Indices::Phase::vapor) << ','
               << water_(state.surfacePhaseRate) << '\n';
    }

    WellControlSwitchOutputOptions options_;
    PetscMPIInt rank_{0};
    mutable std::size_t switchCount_{0};
};

} // namespace MPMC
