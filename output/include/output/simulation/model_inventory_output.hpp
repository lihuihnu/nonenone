/**
 * @file model_inventory_output.hpp
 * @brief 保存各相、溶解、吸附和滞留库存量。
 */
#pragma once

#include <common/console.hpp>

#include <output/core/format.hpp>
#include <output/core/types.hpp>
#include <output/metrics/component_totals.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace MPMC
{

// Forward declaration keeps the output module independent of a concrete grid/runtime header.
template <class Indices>
struct NaturalGlobalInventory;

/** @brief 面向用户的全局质量库存诊断配置。 */
struct ModelInventoryOutputOptions final
{
    std::filesystem::path resultDirectory{"./results"};
    std::filesystem::path historyFile{"model_inventory.csv"};

    bool printState{true};
    bool writeHistory{true};
    bool writeMassTotals{true};
    bool writeInitialState{true};
    std::size_t printEveryOutputSteps{1};
    int dissolvedCO2Component{-1};

    double secondsPerDisplayUnit{86400.0};
    std::string displayTimeUnit{"day"};
};

/**
 * @brief 打印并持久化全局归约后的组分/质量库存。
 *
 * `Runtime::evaluateGlobalInventory()` performs the MPI reduction. This class
 * only owns presentation policy and files, so the physical inventory formulas
 * stay inside the Natural runtime. Rank 0 writes; all ranks may call `write()`
 * with the same arguments.
 */
template <class Indices, class Runtime>
class ModelInventoryOutput final
{
public:
    using Inventory = NaturalGlobalInventory<Indices>;

    /** @brief 在 rank 0 创建 `options` 请求的输出文件。 */
    ModelInventoryOutput(
        ModelInventoryOutputOptions options,
        PetscMPIInt rank)
        : options_(std::move(options)),
          rank_(rank)
    {
        if (options_.printEveryOutputSteps == 0)
            throw std::invalid_argument(
                "ModelInventoryOutput print interval must be positive.");

        if (rank_ != 0)
            return;

        std::filesystem::create_directories(options_.resultDirectory);

        if (options_.writeHistory)
        {
            const auto file = options_.resultDirectory / options_.historyFile;
            std::ofstream stream(file, std::ios::out | std::ios::trunc);
            if (!stream)
                throw std::runtime_error(
                    "Failed to create model inventory history: " + file.string());

            stream << "step,time_" << options_.displayTimeUnit;
            if constexpr (Indices::hasIndependentWaterConservation)
                stream << ",water_mass_kg";
            stream << ",dissolved_co2_mass_kg,trapped_gas_mass_kg";

            for (int c = 0; c < Indices::numComponents; ++c)
                stream << ",fluid_component_" << c << "_kg";
            for (int c = 0; c < Indices::numComponents; ++c)
                stream << ",oil_phase_component_" << c << "_kg";
            for (int c = 0; c < Indices::numComponents; ++c)
                stream << ",trapped_component_" << c << "_kg";
            for (int c = 0; c < Indices::numComponents; ++c)
                stream << ",adsorbed_component_" << c << "_kg";
            stream << '\n';
        }

        if (options_.writeMassTotals)
        {
            for (std::size_t index = 0;
                 index < static_cast<std::size_t>(MassSeries::Count); ++index)
            {
                const auto series = static_cast<MassSeries>(index);
                if (!massSeriesEnabled<Indices>(series))
                    continue;

                std::ofstream stream(
                    options_.resultDirectory / massSeriesInfo(series).totalFilename,
                    std::ios::out | std::ios::trunc);
                writeMassHeader<static_cast<std::size_t>(Indices::numComponents)>(
                    stream, series);
            }
        }
    }

    /**
     * @brief 在一个已接受输出状态计算并写出质量库存。
     *
     * @param step Fixed-output index; zero denotes the initial state.
     * @param timeSeconds Physical time in seconds.
     * @param runtime Natural runtime used for the global inventory reduction.
     * @param solution Accepted primary-variable vector.
     */
    void write(
        std::size_t step,
        double timeSeconds,
        Runtime &runtime,
        Vec solution) const
    {
        const bool initial = step == 0;

        const bool printThisStep =
            options_.printState &&
            ((initial && options_.writeInitialState) ||
             (!initial &&
              step % options_.printEveryOutputSteps == 0));

        const bool writeThisStep =
            (options_.writeHistory ||
             options_.writeMassTotals) &&
            (!initial || options_.writeInitialState);

        if (!printThisStep && !writeThisStep)
            return;

        const Inventory inventory =
            runtime.evaluateGlobalInventory(solution);
        writeInventory(step, timeSeconds, inventory);
    }

    /** @brief 写出已经由 Natural runtime 完成全局归约的库存。 */
    void writeInventory(
        std::size_t step,
        double timeSeconds,
        const Inventory &inventory) const
    {
        const bool initial = step == 0;
        const bool printThisStep =
            options_.printState &&
            ((initial && options_.writeInitialState) ||
             (!initial && step % options_.printEveryOutputSteps == 0));
        const bool writeThisStep =
            (options_.writeHistory || options_.writeMassTotals) &&
            (!initial || options_.writeInitialState);

        if (printThisStep)
            print_(step, timeSeconds, inventory);

        if (writeThisStep)
        {
            if (options_.writeHistory)
                appendCsv_(step, timeSeconds, inventory);

            if (options_.writeMassTotals)
                appendMassTotals_(timeSeconds, inventory);
        }
    }

private:
    static double arrayTotal_(
        const std::array<double, Indices::numComponents> &values)
    {
        return std::accumulate(values.begin(), values.end(), 0.0);
    }

    void print_(
        std::size_t step,
        double timeSeconds,
        const Inventory &inventory) const
    {
        if (rank_ != 0)
            return;

        const std::string header = consoleCentered("MODEL INVENTORY");
        PetscPrintf(
            PETSC_COMM_SELF,
            "\n%s\n"
            "  output step=%zu  time=%.10g %s\n"
            "%s\n",
            header.c_str(),
            step,
            timeSeconds / options_.secondsPerDisplayUnit,
            options_.displayTimeUnit.c_str(),
            consoleRule('-').c_str());

        if constexpr (Indices::hasIndependentWaterConservation)
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "  Water mass                   : %.12e kg\n",
                inventory.waterMass);
        }

        PetscPrintf(
            PETSC_COMM_SELF,
            "  Oil-phase component mass     : %.12e kg\n",
            arrayTotal_(inventory.oilPhaseComponentMass));

        if constexpr (Indices::hasAqueousCO2Dissolution)
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "  Dissolved CO2 mass           : %.12e kg\n",
                inventory.dissolvedCO2Mass);
        }

        if constexpr (Indices::hasLandTrapping)
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "  Land trapped gas mass        : %.12e kg\n",
                inventory.trappedGasMass);
        }

        if constexpr (Indices::hasAdsorption)
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "  Adsorbed component mass      : %.12e kg\n",
                arrayTotal_(inventory.adsorbedComponentMass));
        }

        PetscPrintf(
            PETSC_COMM_SELF,
            "%s\n", consoleRule().c_str());
    }

    void appendCsv_(
        std::size_t step,
        double timeSeconds,
        const Inventory &inventory) const
    {
        if (rank_ != 0)
            return;

        const auto file = options_.resultDirectory / options_.historyFile;
        std::ofstream stream(file, std::ios::out | std::ios::app);
        if (!stream)
            throw std::runtime_error(
                "Failed to open model inventory history: " + file.string());

        stream << std::setprecision(16)
               << step << ','
               << timeSeconds / options_.secondsPerDisplayUnit;
        if constexpr (Indices::hasIndependentWaterConservation)
            stream << ',' << inventory.waterMass;
        stream << ',' << inventory.dissolvedCO2Mass
               << ',' << inventory.trappedGasMass;

        for (double value : inventory.fluidComponentMass)
            stream << ',' << value;
        for (double value : inventory.oilPhaseComponentMass)
            stream << ',' << value;
        for (double value : inventory.trappedComponentMass)
            stream << ',' << value;
        for (double value : inventory.adsorbedComponentMass)
            stream << ',' << value;
        stream << '\n';
    }

    void appendMassTotals_(
        double timeSeconds,
        const Inventory &inventory) const
    {
        if (rank_ != 0)
            return;

        using Totals = ComponentTotals<
            static_cast<std::size_t>(Indices::numComponents)>;

        const auto writeComponents = [&](
            MassSeries series,
            const std::array<double, Indices::numComponents> &componentMass)
        {
            Totals totals;
            totals.component = componentMass;
            std::ofstream stream(
                options_.resultDirectory / massSeriesInfo(series).totalFilename,
                std::ios::out | std::ios::app);
            writeMassTotalLine(stream, series, timeSeconds, totals);
        };

        if constexpr (Indices::hasLandTrapping)
            writeComponents(MassSeries::Trapped, inventory.trappedComponentMass);

        if constexpr (Indices::hasAqueousCO2Dissolution)
        {
            Totals totals;
            if (options_.dissolvedCO2Component >= 0 &&
                options_.dissolvedCO2Component < Indices::numComponents)
            {
                totals.component[static_cast<std::size_t>(options_.dissolvedCO2Component)] =
                    inventory.dissolvedCO2Mass;
            }
            std::ofstream stream(
                options_.resultDirectory /
                    massSeriesInfo(MassSeries::Dissolved).totalFilename,
                std::ios::out | std::ios::app);
            writeMassTotalLine(
                stream, MassSeries::Dissolved, timeSeconds, totals);
        }

        if constexpr (Indices::hasAdsorption)
            writeComponents(MassSeries::Adsorbed, inventory.adsorbedComponentMass);

        writeComponents(
            MassSeries::OilPhaseComponents, inventory.oilPhaseComponentMass);
    }

    ModelInventoryOutputOptions options_;
    PetscMPIInt rank_{0};
};

} // namespace MPMC
