/**
 * @file detailed_well_output.hpp
 * @brief 收敛时间步的井 BHP、流量和穿孔详细诊断。
 */
#pragma once

#include <common/console.hpp>

#include <well/state.hpp>
#include <well/types.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

namespace MPMC
{

struct DetailedWellOutputOptions final
{
    std::filesystem::path resultDirectory{"./results"};
    std::filesystem::path historyFile{"well_history.csv"};

    bool printState{true};
    bool printPhaseDetails{true};
    bool writeHistory{true};
    std::size_t printEveryOutputSteps{1};
    bool writeInitialState{true};

    double secondsPerDisplayUnit{86400.0};
    std::string displayTimeUnit{"day"};
    double pressureScale{1.0e5};
    std::string pressureUnit{"bar"};
    std::vector<std::string> componentNames{};
};

/**
 * @brief 收敛井状态的详细诊断输出。
 *
 * Screen output is intentionally compact at the first level and optionally
 * expands each well with surface/reservoir/mass rates, flow-weighted reservoir
 * density, control limits and total well index.  CSV history is written at
 * every fixed output time so it can be used directly for plots/statistics.
 */
template <class Indices, class Runtime>
class DetailedWellOutput final
{
public:
    using WellStateType = WellState<Indices>;

    DetailedWellOutput(DetailedWellOutputOptions options, PetscMPIInt rank)
        : options_(std::move(options)), rank_(rank)
    {
        if (options_.printEveryOutputSteps == 0)
            throw std::invalid_argument("DetailedWellOutput printEveryOutputSteps must be positive.");

        if (options_.componentNames.empty())
        {
            options_.componentNames.reserve(
                static_cast<std::size_t>(Indices::numComponents));
            for (int component = 0; component < Indices::numComponents; ++component)
                options_.componentNames.emplace_back(
                    "component_" + std::to_string(component));
        }
        if (options_.componentNames.size() !=
            static_cast<std::size_t>(Indices::numComponents))
        {
            throw std::invalid_argument(
                "DetailedWellOutput component-name count does not match Indices.");
        }

        if (rank_ == 0)
            std::filesystem::create_directories(options_.resultDirectory);

        if (rank_ == 0 && options_.writeHistory)
        {
            const auto file = options_.resultDirectory / options_.historyFile;
            std::ofstream stream(file, std::ios::out | std::ios::trunc);
            if (!stream)
                throw std::runtime_error("Failed to create well history file: " + file.string());

            stream
                << "step,time,id,name,type,status,primary_control,primary_target,active_control,control_unit,"
                << "control_target,control_actual,control_error,bhp,bhp_min_limit,bhp_max_limit,water_rate_max_limit,"
                << "q_oil_surface,q_gas_surface,q_water_surface,q_total_surface,"
                << "q_oil_reservoir,q_gas_reservoir,q_water_reservoir,q_total_reservoir,"
                << "m_oil,m_gas,m_water,m_total,";
            for (const auto &name : options_.componentNames)
                stream << "m_component_" << csvToken_(name) << ',';
            stream
                << "rho_oil_reservoir,rho_gas_reservoir,rho_water_reservoir,"
                << "bhp_input_cell,perforation_input_cells,perforations,total_well_index\n";
        }
    }

    void write(std::size_t step, double timeSeconds, Runtime &runtime, Vec solution) const
    {
        if constexpr (!Indices::hasWellUnknown)
            return;
        if (runtime.wells().empty())
            return;

        const bool initial = step == 0;
        const bool printThisStep =
            options_.printState &&
            ((initial && options_.writeInitialState) ||
             (!initial && step % options_.printEveryOutputSteps == 0));
        const bool writeThisStep =
            options_.writeHistory && (!initial || options_.writeInitialState);

        if (!printThisStep && !writeThisStep)
            return;

        const auto states = runtime.evaluateWellStates(solution);
        if (printThisStep)
            printStates_(step, timeSeconds, runtime, states);
        if (writeThisStep)
            appendHistory_(step, timeSeconds, runtime, states);
    }

private:
    [[nodiscard]] static std::string csvToken_(std::string name)
    {
        for (char &character : name)
        {
            const auto value = static_cast<unsigned char>(character);
            if (!std::isalnum(value) && character != '_')
                character = '_';
        }
        return name;
    }

    [[nodiscard]] static const char *activityName_(bool active) noexcept
    {
        return active ? "ACTIVE" : "SHUT";
    }

    [[nodiscard]] static double phase_(const std::array<double, Indices::numPhases> &values, int phase)
    {
        return values[static_cast<std::size_t>(phase)];
    }

    template <class Array>
    [[nodiscard]] static double oil_(const Array &values) { return phase_(values, Indices::Phase::liquid); }
    template <class Array>
    [[nodiscard]] static double gas_(const Array &values) { return phase_(values, Indices::Phase::vapor); }
    template <class Array>
    [[nodiscard]] static double water_(const Array &values)
    {
        if constexpr (Indices::hasWater)
            return phase_(values, Indices::Phase::water);
        return 0.0;
    }

    [[nodiscard]] static double total_(const std::array<double, Indices::numPhases> &values)
    {
        double result = 0.0;
        for (double value : values)
            result += value;
        return result;
    }

    template <class Well>
    [[nodiscard]] double scaledTarget_(WellControl control, double target) const
    {
        return control == WellControl::Bhp ? target / options_.pressureScale : target;
    }

    template <class Well>
    [[nodiscard]] double controlTarget_(const Well &well) const
    {
        return scaledTarget_<Well>(well.control, well.target);
    }

    template <class Well>
    [[nodiscard]] double primaryTarget_(const Well &well) const
    {
        return scaledTarget_<Well>(well.primaryControl, well.primaryTarget);
    }

    template <class Well>
    [[nodiscard]] double controlActual_(const Well &well, const WellStateType &state) const
    {
        if (well.control == WellControl::Bhp)
            return state.bottomHolePressure / options_.pressureScale;
        return state.controlledRateMagnitude(well.type, well.control);
    }

    template <class Well>
    [[nodiscard]] double controlError_(const Well &well, const WellStateType &state) const
    {
        return controlActual_(well, state) - controlTarget_(well);
    }

    template <class Well>
    [[nodiscard]] const char *controlUnit_(const Well &well) const noexcept
    {
        return well.control == WellControl::Bhp ? options_.pressureUnit.c_str() : "m3/s";
    }

    template <class Well>
    [[nodiscard]] static double totalWellIndex_(const Well &well)
    {
        double result = 0.0;
        for (const auto &perf : well.perforations)
            result += perf.wellIndex;
        return result;
    }

    template <class Well>
    [[nodiscard]] double minBhp_(const Well &well) const
    {
        return well.limits.minimumBhp
            ? *well.limits.minimumBhp / options_.pressureScale
            : std::numeric_limits<double>::quiet_NaN();
    }

    template <class Well>
    [[nodiscard]] double maxBhp_(const Well &well) const
    {
        return well.limits.maximumBhp
            ? *well.limits.maximumBhp / options_.pressureScale
            : std::numeric_limits<double>::quiet_NaN();
    }

    template <class Well>
    [[nodiscard]] static double maxWaterRate_(const Well &well)
    {
        return well.limits.maximumWaterRate
            ? *well.limits.maximumWaterRate
            : std::numeric_limits<double>::quiet_NaN();
    }

    template <class Well>
    [[nodiscard]] static std::string inputCellList_(
        const Runtime &runtime,
        const Well &well)
    {
        std::ostringstream stream;
        for (std::size_t i = 0; i < well.perforations.size(); ++i)
        {
            if (i != 0)
                stream << ';';
            stream << runtime.inputCellIndex(
                well.perforations[i].currentCellId);
        }
        return stream.str();
    }

    void printStates_(
        std::size_t step,
        double timeSeconds,
        const Runtime &runtime,
        const std::vector<WellStateType> &states) const
    {
        if (rank_ != 0 || !options_.printState)
            return;

        const auto &wells = runtime.wells();
        if (wells.size() != states.size())
            throw std::runtime_error("Well state count does not match runtime wells.");

        const std::string header = consoleCentered("WELL STATE");
        const std::string separator = consoleRule('-');
        PetscPrintf(
            PETSC_COMM_SELF,
            "\n%s\n"
            "  output step=%zu  time=%.10g %s  sign: injection(+), production(-)\n"
            "%s\n"
            "%-4s %-12s %-9s %-7s %-11s %12s %12s %12s %12s %5s\n",
            header.c_str(),
            step,
            timeSeconds / options_.secondsPerDisplayUnit,
            options_.displayTimeUnit.c_str(),
            separator.c_str(),
            "ID", "NAME", "TYPE", "STATUS", "CONTROL",
            "TARGET", "ACTUAL", "ERROR", "BHP(bar)", "NPERF");

        for (std::size_t i = 0; i < wells.size(); ++i)
        {
            const auto &well = wells[i];
            const auto &state = states[i];
            PetscPrintf(
                PETSC_COMM_SELF,
                "%-4d %-12s %-9s %-7s %-11s %12.5g %12.5g %12.5g %12.5g %5zu\n",
                well.id,
                well.name.c_str(),
                std::string(wellTypeName(well.type)).c_str(),
                activityName_(well.isActive(timeSeconds)),
                std::string(wellControlName(well.control)).c_str(),
                controlTarget_(well),
                controlActual_(well, state),
                controlError_(well, state),
                state.bottomHolePressure / options_.pressureScale,
                well.perforations.size());

            if (options_.printPhaseDetails)
            {
                PetscPrintf(
                    PETSC_COMM_SELF,
                    "      primary=%s target=%.6g  active-unit=%s  WI(total)=%.6e\n"
                    "      input-cells: bhp=%" PetscInt_FMT "  perf=[%s]\n"
                    "      limits: BHP[min/max]=%.6g / %.6g bar, water-rate-max=%.6g m3/s\n"
                    "      q_surface   [O/G/W/T] = % .6e  % .6e  % .6e  % .6e m3/s\n"
                    "      q_reservoir [O/G/W/T] = % .6e  % .6e  % .6e  % .6e m3/s\n"
                    "      mass rate   [O/G/W/T] = % .6e  % .6e  % .6e  % .6e kg/s\n"
                    "      rho_res     [O/G/W]   = % .6e  % .6e  % .6e kg/m3\n",
                    std::string(wellControlName(well.primaryControl)).c_str(),
                    primaryTarget_(well),
                    controlUnit_(well),
                    totalWellIndex_(well),
                    runtime.inputCellIndex(well.bhpCellId),
                    inputCellList_(runtime, well).c_str(),
                    minBhp_(well), maxBhp_(well), maxWaterRate_(well),
                    oil_(state.surfacePhaseRate), gas_(state.surfacePhaseRate), water_(state.surfacePhaseRate), total_(state.surfacePhaseRate),
                    oil_(state.reservoirPhaseRate), gas_(state.reservoirPhaseRate), water_(state.reservoirPhaseRate), total_(state.reservoirPhaseRate),
                    oil_(state.phaseMassRate), gas_(state.phaseMassRate), water_(state.phaseMassRate), total_(state.phaseMassRate),
                    oil_(state.flowWeightedDensity), gas_(state.flowWeightedDensity), water_(state.flowWeightedDensity));
                std::ostringstream components;
                components << "      component mass rate [kg/s] =";
                for (std::size_t component = 0;
                     component < state.componentMassRate.size(); ++component)
                {
                    components << ' ' << options_.componentNames[component]
                               << '=' << std::scientific
                               << std::setprecision(6)
                               << state.componentMassRate[component];
                }
                components << '\n';
                PetscPrintf(PETSC_COMM_SELF, "%s", components.str().c_str());
            }
        }

        PetscPrintf(
            PETSC_COMM_SELF,
            "%s\n",
            consoleRule().c_str());
    }

    void appendHistory_(
        std::size_t step,
        double timeSeconds,
        const Runtime &runtime,
        const std::vector<WellStateType> &states) const
    {
        if (rank_ != 0 || !options_.writeHistory)
            return;

        const auto &wells = runtime.wells();
        if (wells.size() != states.size())
            throw std::runtime_error("Well state count does not match runtime wells.");

        const auto file = options_.resultDirectory / options_.historyFile;
        std::ofstream stream(file, std::ios::out | std::ios::app);
        if (!stream)
            throw std::runtime_error("Failed to open well history file: " + file.string());

        stream << std::setprecision(16);
        const double displayTime = timeSeconds / options_.secondsPerDisplayUnit;

        for (std::size_t i = 0; i < wells.size(); ++i)
        {
            const auto &well = wells[i];
            const auto &state = states[i];
            stream
                << step << ',' << displayTime << ','
                << well.id << ',' << well.name << ',' << wellTypeName(well.type) << ','
                << activityName_(well.isActive(timeSeconds)) << ','
                << wellControlName(well.primaryControl) << ',' << primaryTarget_(well) << ','
                << wellControlName(well.control) << ',' << controlUnit_(well) << ','
                << controlTarget_(well) << ',' << controlActual_(well, state) << ',' << controlError_(well, state) << ','
                << state.bottomHolePressure / options_.pressureScale << ','
                << minBhp_(well) << ',' << maxBhp_(well) << ',' << maxWaterRate_(well) << ','
                << oil_(state.surfacePhaseRate) << ',' << gas_(state.surfacePhaseRate) << ',' << water_(state.surfacePhaseRate) << ',' << total_(state.surfacePhaseRate) << ','
                << oil_(state.reservoirPhaseRate) << ',' << gas_(state.reservoirPhaseRate) << ',' << water_(state.reservoirPhaseRate) << ',' << total_(state.reservoirPhaseRate) << ','
                << oil_(state.phaseMassRate) << ',' << gas_(state.phaseMassRate) << ',' << water_(state.phaseMassRate) << ',' << total_(state.phaseMassRate) << ',';
            for (double rate : state.componentMassRate)
                stream << rate << ',';
            stream
                << oil_(state.flowWeightedDensity) << ',' << gas_(state.flowWeightedDensity) << ',' << water_(state.flowWeightedDensity) << ','
                << runtime.inputCellIndex(well.bhpCellId) << ','
                << inputCellList_(runtime, well) << ','
                << well.perforations.size() << ',' << totalWellIndex_(well) << '\n';
        }
    }

    DetailedWellOutputOptions options_;
    PetscMPIInt rank_{0};
};

} // namespace MPMC
