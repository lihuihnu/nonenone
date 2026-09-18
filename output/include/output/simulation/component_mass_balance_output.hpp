/**
 * @file component_mass_balance_output.hpp
 * @brief 逐时间步保存逐组分质量守恒诊断。
 */
#pragma once

#include <common/console.hpp>

#include <output/metrics/component_mass_balance.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief 离散全局组分守恒诊断的展示策略。
 *
 * The diagnostic assumes the reservoir's external boundaries are closed.  All
 * current MPMC_SCW standard reservoir cases use no-flow external boundaries;
 * if open boundary conditions are introduced later, their component fluxes
 * must be added to the ledger in the same accepted-step hook as well sources.
 */
struct ComponentMassBalanceOutputOptions final
{
    std::filesystem::path resultDirectory{"./results"};
    std::filesystem::path historyFile{"component_mass_balance.csv"};

    bool enabled{true};
    bool printState{true};
    bool writeHistory{true};
    bool writeInitialState{true};
    bool auditInternalFaceConservation{false};
    std::size_t printEveryOutputSteps{1};

    std::vector<std::string> componentNames{};
    double secondsPerDisplayUnit{86400.0};
    std::string displayTimeUnit{"day"};
};

/**
 * @brief 离散 MPI 全局逐组分守恒诊断器。
 *
 * Reservoir inventory is evaluated as
 *
 *   M_i = M_i,fluid + M_i,adsorbed,
 *
 * because adsorption is an additional storage term in the component equation.
 * For the legacy formulation, the independent H2O conservation equation is
 * appended as one additional conserved-species row using inventory.waterMass.
 * Land-trapped gas is already part of the gas-phase fluid inventory and is
 * deliberately not added again.
 *
 * The accepted-step hook accumulates the exact well source used by the fully
 * implicit residual at the accepted end-of-step state.  Thus the diagnostic
 * checks the same backward-Euler discrete balance solved by SNES rather than a
 * post-processed approximation based only on fixed output times.
 */
template <class Indices, class Runtime>
class ComponentMassBalanceOutput final
{
public:
    static constexpr std::size_t componentCount =
        static_cast<std::size_t>(Indices::numComponents);
    static constexpr std::size_t N =
        componentCount +
        (Indices::hasIndependentWaterConservation ? std::size_t{1} : std::size_t{0});
    using Array = std::array<double, N>;
    using Ledger = ComponentMassBalanceLedger<N>;
    using Snapshot = ComponentMassBalanceSnapshot<N>;

    ComponentMassBalanceOutput(
        ComponentMassBalanceOutputOptions options,
        PetscMPIInt rank)
        : options_(std::move(options)),
          rank_(rank)
    {
        if (!options_.enabled)
            return;
        if (options_.printEveryOutputSteps == 0)
            throw std::invalid_argument(
                "ComponentMassBalanceOutput print interval must be positive.");
        if (options_.secondsPerDisplayUnit <= 0.0 ||
            !std::isfinite(options_.secondsPerDisplayUnit))
            throw std::invalid_argument(
                "ComponentMassBalanceOutput display-time scale must be positive.");

        normalizeComponentNames_();

        if (rank_ != 0)
            return;

        std::filesystem::create_directories(options_.resultDirectory);
        if (options_.writeHistory)
        {
            const auto file = options_.resultDirectory / options_.historyFile;
            std::ofstream stream(file, std::ios::out | std::ios::trunc);
            if (!stream)
                throw std::runtime_error(
                    "Failed to create component mass-balance history: " +
                    file.string());

            stream
                << "step,time_" << options_.displayTimeUnit
                << ",component_index,component"
                << ",initial_inventory_kg,current_inventory_kg"
                << ",cumulative_injected_kg,cumulative_produced_kg"
                << ",expected_inventory_kg,balance_error_kg,relative_error";
            if (options_.auditInternalFaceConservation)
            {
                stream
                    << ",internal_face_flux_imbalance_kg_s"
                    << ",cumulative_internal_face_imbalance_kg";
            }
            stream << '\n';
        }
    }

    /**
     * @brief 记录一个已接受的内部时间步。
     *
     * This must be called for every accepted adaptive/fixed internal step, not
     * only fixed output states.  Rejected attempts intentionally never call it.
     */
    void acceptedStep(Runtime &runtime, Vec solution, double acceptedTimeSeconds)
    {
        if (!options_.enabled)
            return;
        if (!ledger_.initialized())
            throw std::logic_error(
                "Component mass-balance initial state was not written before the first accepted step.");

        const auto rates = runtime.evaluateGlobalWellComponentRates(solution);
        Array injected{};
        Array produced{};
        for (std::size_t c = 0; c < componentCount; ++c)
        {
            injected[c] = rates.injected[c];
            produced[c] = rates.produced[c];
        }
        if constexpr (Indices::hasIndependentWaterConservation)
        {
            injected[componentCount] = rates.waterInjected;
            produced[componentCount] = rates.waterProduced;
        }

        const double dt =
            acceptedTimeSeconds - ledger_.lastAcceptedTime();

        if (options_.auditInternalFaceConservation)
        {
            const auto face =
                runtime.evaluateGlobalInternalFaceFluxImbalance(solution);
            for (std::size_t comp = 0; comp < componentCount; ++comp)
            {
                lastInternalFaceImbalance_[comp] =
                    face.component[comp];
                cumulativeInternalFaceImbalance_[comp] +=
                    face.component[comp] * dt;
            }
            if constexpr (Indices::hasIndependentWaterConservation)
            {
                lastInternalFaceImbalance_[componentCount] =
                    face.independentWater;
                cumulativeInternalFaceImbalance_[componentCount] +=
                    face.independentWater * dt;
            }
        }

        ledger_.accept(acceptedTimeSeconds, injected, produced);
    }

    /** @brief 在一个已接受固定输出状态打印/写出守恒误差。 */
    void write(
        std::size_t step,
        double timeSeconds,
        Runtime &runtime,
        Vec solution)
    {
        if (!options_.enabled)
            return;

        const auto inventory = runtime.evaluateGlobalInventory(solution);
        writeInventory(step, timeSeconds, inventory);
    }

    /** @brief 使用已完成一次 MPI 归约的库存打印/写出守恒误差。 */
    template <class Inventory>
    void writeInventory(
        std::size_t step,
        double timeSeconds,
        const Inventory &inventory)
    {
        if (!options_.enabled)
            return;

        const bool initial = step == 0;
        const Array current = conservedInventory_(inventory);

        if (initial)
        {
            if (ledger_.initialized())
                throw std::logic_error(
                    "Component mass-balance initial state was initialized more than once.");
            ledger_.initialize(current, timeSeconds);
        }
        else
        {
            if (!ledger_.initialized())
                throw std::logic_error(
                    "Component mass-balance output requested before initialization.");

            const double tolerance =
                1.0e-10 * std::max(1.0, std::abs(timeSeconds));
            if (std::abs(ledger_.lastAcceptedTime() - timeSeconds) > tolerance)
                throw std::logic_error(
                    "Component mass-balance ledger time does not match the accepted output time.");
        }

        const Snapshot snapshot = ledger_.snapshot(current);
        const bool printThisStep =
            options_.printState &&
            ((initial && options_.writeInitialState) ||
             (!initial && step % options_.printEveryOutputSteps == 0));
        const bool writeThisStep =
            options_.writeHistory && (!initial || options_.writeInitialState);

        if (printThisStep)
            print_(step, timeSeconds, snapshot);
        if (writeThisStep)
            appendCsv_(step, timeSeconds, snapshot);
    }

private:
    void normalizeComponentNames_()
    {
        if (options_.componentNames.empty())
        {
            options_.componentNames.reserve(N);
            for (std::size_t c = 0; c < N; ++c)
                options_.componentNames.emplace_back("component_" + std::to_string(c));
            return;
        }
        if constexpr (Indices::hasIndependentWaterConservation)
        {
            // Case configs naturally list compositional EOS components only.
            // Append the independently conserved H2O equation automatically.
            if (options_.componentNames.size() == componentCount)
                options_.componentNames.emplace_back("H2O");
        }
        if (options_.componentNames.size() != N)
            throw std::invalid_argument(
                "ComponentMassBalanceOutput conserved-species name count does not match Indices.");
    }

    template <class Inventory>
    static Array conservedInventory_(const Inventory &inventory)
    {
        Array result{};
        for (std::size_t c = 0; c < componentCount; ++c)
        {
            // Adsorption is a genuine storage term. Land trapped mass is only a
            // subset of fluid gas mass and therefore must not be added here.
            result[c] = inventory.fluidComponentMass[c] +
                        inventory.adsorbedComponentMass[c];
        }
        if constexpr (Indices::hasIndependentWaterConservation)
            result[componentCount] = inventory.waterMass;
        return result;
    }

    static double sum_(const Array &values)
    {
        return std::accumulate(values.begin(), values.end(), 0.0);
    }

    static double maxAbs_(const Array &values)
    {
        double result = 0.0;
        for (double value : values)
            result = std::max(result, std::abs(value));
        return result;
    }

    void print_(
        std::size_t step,
        double timeSeconds,
        const Snapshot &s) const
    {
        if (rank_ != 0)
            return;

        const std::string header = consoleCentered("COMPONENT MASS BALANCE");
        PetscPrintf(
            PETSC_COMM_SELF,
            "\n%s\n"
            "  output step=%zu  time=%.10g %s\n"
            "  closed external boundaries; cumulative well terms use accepted BE steps\n"
            "%s\n"
            "%-10s %14s %14s %12s %12s %12s %11s\n",
            header.c_str(),
            step,
            timeSeconds / options_.secondsPerDisplayUnit,
            options_.displayTimeUnit.c_str(),
            consoleRule('-').c_str(),
            "Component", "M0(kg)", "M(kg)", "CumIn(kg)", "CumOut(kg)",
            "Error(kg)", "RelError");

        for (std::size_t c = 0; c < N; ++c)
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "%-10s %14.6e %14.6e %12.4e %12.4e %12.4e %11.3e\n",
                options_.componentNames[c].c_str(),
                s.initial[c],
                s.current[c],
                s.cumulativeInjected[c],
                s.cumulativeProduced[c],
                s.error[c],
                s.relativeError[c]);
        }

        const double initialTotal = sum_(s.initial);
        const double currentTotal = sum_(s.current);
        const double injectedTotal = sum_(s.cumulativeInjected);
        const double producedTotal = sum_(s.cumulativeProduced);
        const double errorTotal = sum_(s.error);
        const double totalScale = std::max({
            1.0,
            std::abs(initialTotal),
            std::abs(currentTotal),
            injectedTotal + producedTotal});

        PetscPrintf(
            PETSC_COMM_SELF,
            "%-10s %14.6e %14.6e %12.4e %12.4e %12.4e %11.3e\n"
            "  Max |component relative error| = %.6e\n"
            "%s\n",
            "TOTAL",
            initialTotal,
            currentTotal,
            injectedTotal,
            producedTotal,
            errorTotal,
            errorTotal / totalScale,
            maxAbs_(s.relativeError),
            consoleRule().c_str());
    }

    void appendCsv_(
        std::size_t step,
        double timeSeconds,
        const Snapshot &s) const
    {
        if (rank_ != 0)
            return;

        const auto file = options_.resultDirectory / options_.historyFile;
        std::ofstream stream(file, std::ios::out | std::ios::app);
        if (!stream)
            throw std::runtime_error(
                "Failed to open component mass-balance history: " +
                file.string());

        stream << std::setprecision(16);
        for (std::size_t c = 0; c < N; ++c)
        {
            stream
                << step << ','
                << timeSeconds / options_.secondsPerDisplayUnit << ','
                << c << ',' << options_.componentNames[c] << ','
                << s.initial[c] << ','
                << s.current[c] << ','
                << s.cumulativeInjected[c] << ','
                << s.cumulativeProduced[c] << ','
                << s.expected[c] << ','
                << s.error[c] << ','
                << s.relativeError[c];
            if (options_.auditInternalFaceConservation)
            {
                stream
                    << ',' << lastInternalFaceImbalance_[c]
                    << ',' << cumulativeInternalFaceImbalance_[c];
            }
            stream << '\n';
        }
    }

    ComponentMassBalanceOutputOptions options_;
    PetscMPIInt rank_{0};
    Ledger ledger_{};
    Array lastInternalFaceImbalance_{};
    Array cumulativeInternalFaceImbalance_{};
};

} // namespace MPMC
