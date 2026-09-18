/**
 * @file producer_composition_output.hpp
 * @brief 生产井一级组分输出：瞬时质量分数、累计产量、RF 与轻重富集。
 */
#pragma once

#include <output/metrics/producer_composition.hpp>

#include <well/state.hpp>
#include <well/types.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
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
#include <utility>
#include <vector>

namespace MPMC
{

struct ProducerCompositionOutputOptions final
{
    std::filesystem::path resultDirectory{"./results"};
    std::filesystem::path historyFile{"producer_composition.csv"};

    bool enabled{false};
    bool writeHistory{true};
    bool writeInitialState{true};

    std::vector<std::string> componentNames{};
    std::vector<int> lightComponents{};
    std::vector<int> heavyComponents{};
    int waterComponent{-1};

    // Formal slab runs must supply the measured operating-condition PV.
    // Legacy/numerical cases may leave it NaN; cumulative injected volume is
    // still written, but PVI remains NaN rather than using a hidden PV.
    double effectivePoreVolumeM3{
        std::numeric_limits<double>::quiet_NaN()};

    double secondsPerDisplayUnit{1.0};
    std::string displayTimeUnit{"s"};
};

/**
 * @brief 按 accepted internal step 积分每口生产井的组分采出历史。
 *
 * Instantaneous composition is derived from the exact conserved-component
 * well source.  Cumulative production is integrated with accepted end-step
 * rates, matching the Backward-Euler mass-balance convention.
 */
template <class Indices, class Runtime>
class ProducerCompositionOutput final
{
public:
    static constexpr std::size_t N =
        static_cast<std::size_t>(Indices::numComponents);
    using Array = std::array<double, N>;
    using Ledger = ProducerCompositionLedger<N>;
    using WellStateType = WellState<Indices>;

    ProducerCompositionOutput(
        ProducerCompositionOutputOptions options,
        PetscMPIInt rank)
        : options_(std::move(options)), rank_(rank)
    {
        if (!options_.enabled)
            return;

        normalizeComponentNames_();
        validateGroups_();

        if (std::isfinite(options_.effectivePoreVolumeM3) &&
            !(options_.effectivePoreVolumeM3 > 0.0))
        {
            throw std::invalid_argument(
                "Producer-composition effective pore volume must be positive when supplied.");
        }

        if (rank_ != 0)
            return;

        std::filesystem::create_directories(options_.resultDirectory);
        if (options_.writeHistory)
        {
            const auto file =
                options_.resultDirectory / options_.historyFile;
            std::ofstream out(file, std::ios::out | std::ios::trunc);
            if (!out)
                throw std::runtime_error(
                    "Failed to create producer composition history.");

            out << "step,time_" << options_.displayTimeUnit
                << ",producer_id,producer_name"
                << ",cumulative_injected_reservoir_m3,pvi";
            for (const auto &name : options_.componentNames)
            {
                const auto token = csvToken_(name);
                out << ",m_dot_" << token << "_produced_kg_s"
                    << ",Y_" << token << "_mass_fraction"
                    << ",cum_" << token << "_produced_kg"
                    << ",RF_" << token;
            }
            out << ",RF_total_hydrocarbon"
                << ",E_L_over_H_instant_mass"
                << ",E_L_over_H_cumulative_mass\n";
        }
    }

    void write(
        std::size_t step,
        double timeSeconds,
        Runtime &runtime,
        Vec solution)
    {
        if (!options_.enabled)
            return;

        if (!initialized_)
            initialize_(runtime, solution, timeSeconds);

        if (step == 0 && !options_.writeInitialState)
            return;
        if (!options_.writeHistory || rank_ != 0)
            return;

        const auto states = runtime.evaluateWellStates(solution);
        const auto &wells = runtime.wells();
        if (states.size() != wells.size())
            throw std::runtime_error(
                "Producer composition well-state count mismatch.");

        const double pvi =
            std::isfinite(options_.effectivePoreVolumeM3)
                ? cumulativeInjectedReservoirM3_ /
                      options_.effectivePoreVolumeM3
                : std::numeric_limits<double>::quiet_NaN();

        const auto file =
            options_.resultDirectory / options_.historyFile;
        std::ofstream out(file, std::ios::out | std::ios::app);
        if (!out)
            throw std::runtime_error(
                "Failed to open producer composition history.");
        out << std::setprecision(17);

        for (std::size_t i = 0; i < wells.size(); ++i)
        {
            if (wells[i].type != WellType::Producer)
                continue;

            const Array rate = productionMagnitude_(states[i]);
            const auto snapshot = ledgers_[i].snapshot(
                rate, options_.lightComponents, options_.heavyComponents);

            out << step << ','
                << timeSeconds / options_.secondsPerDisplayUnit << ','
                << wells[i].id << ',' << wells[i].name << ','
                << cumulativeInjectedReservoirM3_ << ',' << pvi;

            for (std::size_t c = 0; c < N; ++c)
            {
                out << ',' << snapshot.instantaneousProducedRateKgPerS[c]
                    << ',' << snapshot.instantaneousMassFraction[c]
                    << ',' << snapshot.cumulativeProducedKg[c]
                    << ',' << snapshot.recoveryFraction[c];
            }
            out << ',' << totalHydrocarbonRecovery_(snapshot)
                << ',' << snapshot.instantaneousLightHeavyEnrichment
                << ',' << snapshot.cumulativeLightHeavyEnrichment
                << '\n';
        }
    }

    void acceptedStep(
        Runtime &runtime,
        Vec solution,
        double acceptedTimeSeconds)
    {
        if (!options_.enabled)
            return;
        if (!initialized_)
            throw std::logic_error(
                "Producer composition output must be initialized at step 0.");

        const auto states = runtime.evaluateWellStates(solution);
        const auto &wells = runtime.wells();
        if (states.size() != wells.size())
            throw std::runtime_error(
                "Producer composition well-state count mismatch.");

        const double dt = acceptedTimeSeconds - lastAcceptedTimeSeconds_;
        if (!(dt > 0.0))
            throw std::invalid_argument(
                "Producer composition accepted times must increase.");

        double injectedReservoirRate = 0.0;
        for (std::size_t i = 0; i < wells.size(); ++i)
        {
            if (!wells[i].isActive(acceptedTimeSeconds))
                continue;

            if (wells[i].type == WellType::Injector)
            {
                double rate = 0.0;
                for (double value : states[i].reservoirPhaseRate)
                    rate += value;
                injectedReservoirRate += std::max(0.0, rate);
                continue;
            }

            if (wells[i].type == WellType::Producer)
                ledgers_[i].accept(
                    acceptedTimeSeconds,
                    productionMagnitude_(states[i]));
        }

        cumulativeInjectedReservoirM3_ +=
            injectedReservoirRate * dt;
        lastAcceptedTimeSeconds_ = acceptedTimeSeconds;
    }

private:
    static std::string csvToken_(std::string name)
    {
        for (char &ch : name)
        {
            const unsigned char value =
                static_cast<unsigned char>(ch);
            if (!std::isalnum(value) && ch != '_')
                ch = '_';
        }
        return name;
    }

    void normalizeComponentNames_()
    {
        if (options_.componentNames.empty())
        {
            options_.componentNames.reserve(N);
            for (std::size_t c = 0; c < N; ++c)
                options_.componentNames.emplace_back(
                    "component_" + std::to_string(c));
        }
        if (options_.componentNames.size() != N)
            throw std::invalid_argument(
                "Producer-composition component-name count mismatch.");
    }

    [[nodiscard]] double totalHydrocarbonRecovery_(
        const ProducerCompositionSnapshot<N> &snapshot) const
    {
        double initial = 0.0;
        double produced = 0.0;
        for (std::size_t component = 0; component < N; ++component)
        {
            if (static_cast<int>(component) == options_.waterComponent)
                continue;
            initial += snapshot.initialInventoryKg[component];
            produced += snapshot.cumulativeProducedKg[component];
        }
        return initial > 0.0
            ? produced / initial
            : std::numeric_limits<double>::quiet_NaN();
    }

    void validateGroups_() const
    {
        const auto validate = [](const std::vector<int> &group) {
            for (int component : group)
                if (component < 0 ||
                    component >= static_cast<int>(N))
                    throw std::out_of_range(
                        "Producer-composition group index out of range.");
        };
        validate(options_.lightComponents);
        validate(options_.heavyComponents);
        if (options_.waterComponent < -1 ||
            options_.waterComponent >= static_cast<int>(N))
            throw std::out_of_range(
                "Producer-composition water component index out of range.");
    }

    static Array productionMagnitude_(
        const WellStateType &state)
    {
        Array result{};
        for (std::size_t c = 0; c < N; ++c)
            result[c] =
                std::max(0.0, -state.componentMassRate[c]);
        return result;
    }

    void initialize_(
        Runtime &runtime,
        Vec solution,
        double timeSeconds)
    {
        const auto inventory =
            runtime.evaluateGlobalInventory(solution);
        Array initial{};
        for (std::size_t c = 0; c < N; ++c)
        {
            initial[c] =
                inventory.fluidComponentMass[c] +
                inventory.trappedComponentMass[c] +
                inventory.adsorbedComponentMass[c];
        }

        ledgers_.assign(runtime.wells().size(), Ledger{});
        for (std::size_t i = 0; i < runtime.wells().size(); ++i)
            if (runtime.wells()[i].type == WellType::Producer)
                ledgers_[i].initialize(initial, timeSeconds);

        lastAcceptedTimeSeconds_ = timeSeconds;
        cumulativeInjectedReservoirM3_ = 0.0;
        initialized_ = true;
    }

    ProducerCompositionOutputOptions options_;
    PetscMPIInt rank_{0};
    std::vector<Ledger> ledgers_{};
    double lastAcceptedTimeSeconds_{0.0};
    double cumulativeInjectedReservoirM3_{0.0};
    bool initialized_{false};
};

} // namespace MPMC
