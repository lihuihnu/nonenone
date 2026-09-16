/**
 * @file reservoir_diagnostics.hpp
 * @brief 汇总压力、饱和度、相态和储层物性诊断。
 */
#pragma once

#include <common/console.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>

namespace MPMC
{

struct ReservoirDiagnosticsOutputOptions final
{
    std::filesystem::path resultDirectory{"./results"};
    std::filesystem::path historyFile{"reservoir_diagnostics.csv"};
    bool printState{true};
    bool writeHistory{true};
    std::size_t printEveryOutputSteps{1};
    bool writeInitialState{true};
    double secondsPerDisplayUnit{86400.0};
    std::string displayTimeUnit{"day"};
    double pressureScale{1.0e5};
    std::string pressureUnit{"bar"};
};

/**
 * @brief Natural runtime 的全局储层状态诊断输出。
 *
 * The Runtime only needs to provide evaluateGlobalDiagnostics(solution).
 * All quantities are MPI-global and therefore safe to print/write on rank 0.
 */
template <class Indices, class Runtime>
class ReservoirDiagnosticsOutput final
{
public:
    ReservoirDiagnosticsOutput(
        ReservoirDiagnosticsOutputOptions options,
        PetscMPIInt rank)
        : options_(std::move(options)), rank_(rank)
    {
        if (options_.printEveryOutputSteps == 0)
            throw std::invalid_argument("Reservoir diagnostics print frequency must be positive.");

        if (rank_ == 0)
            std::filesystem::create_directories(options_.resultDirectory);

        if (rank_ == 0 && options_.writeHistory)
        {
            const auto file = options_.resultDirectory / options_.historyFile;
            std::ofstream stream(file, std::ios::out | std::ios::trunc);
            if (!stream)
                throw std::runtime_error("Failed to create reservoir diagnostics file: " + file.string());

            stream
                << "step,time,cells,bulk_volume,pore_volume,"
                << "pressure_min,pressure_avg,pressure_max,"
                << "porosity_min,porosity_avg,porosity_max,"
                << "s_o_min,s_o_avg,s_o_max,s_g_min,s_g_avg,s_g_max,s_w_min,s_w_avg,s_w_max,"
                << "rho_o_min,rho_o_avg,rho_o_max,rho_g_min,rho_g_avg,rho_g_max,rho_w_min,rho_w_avg,rho_w_max,"
                << "mu_o_min,mu_o_avg,mu_o_max,mu_g_min,mu_g_avg,mu_g_max,mu_w_min,mu_w_avg,mu_w_max,";

            if constexpr (Indices::fullyCompositionalThreePhase)
            {
                stream << "oil_only_cells,gas_only_cells,water_only_cells,"
                       << "oil_gas_cells,oil_water_cells,gas_water_cells,oil_gas_water_cells,"
                       << "suppressed_oil_cells,suppressed_gas_cells,suppressed_water_cells,";
            }
            else
            {
                stream << "liquid_only_cells,vapor_only_cells,two_phase_cells,";
            }

            stream
                << "aqueous_co2_mass_fraction_avg,aqueous_co2_mass_fraction_max,"
                << "trapped_gas_saturation_avg,trapped_gas_saturation_max\n";
        }
    }

    void write(std::size_t step, double timeSeconds, Runtime &runtime, Vec solution) const
    {
        const bool initial = step == 0;
        const bool printThisStep =
            options_.printState &&
            ((initial && options_.writeInitialState) ||
             (!initial && step % options_.printEveryOutputSteps == 0));
        const bool writeThisStep =
            options_.writeHistory && (!initial || options_.writeInitialState);

        if (!printThisStep && !writeThisStep)
            return;

        const auto diagnostics = runtime.evaluateGlobalDiagnostics(solution);
        if (printThisStep)
            print_(step, timeSeconds, diagnostics);
        if (writeThisStep)
            append_(step, timeSeconds, diagnostics);
    }

private:
    template <class Diagnostics>
    void print_(std::size_t step, double timeSeconds, const Diagnostics &d) const
    {
        if (rank_ != 0)
            return;

        const std::string header = consoleCentered("RESERVOIR DIAGNOSTICS");
        PetscPrintf(
            PETSC_COMM_SELF,
            "\n%s\n"
            "  output step=%zu  time=%.10g %s\n"
            "%s\n"
            "  Cells / bulk volume / pore volume : %lld / %.6e / %.6e m3\n"
            "  Pressure [min/avg/max]             : %.6g / %.6g / %.6g %s\n"
            "  Porosity [min/avg/max]             : %.6g / %.6g / %.6g\n",
            header.c_str(),
            step,
            timeSeconds / options_.secondsPerDisplayUnit,
            options_.displayTimeUnit.c_str(),
            consoleRule('-').c_str(),
            d.cellCount,
            d.bulkVolume,
            d.poreVolume,
            d.pressureMinimum / options_.pressureScale,
            d.pressureAverage / options_.pressureScale,
            d.pressureMaximum / options_.pressureScale,
            options_.pressureUnit.c_str(),
            d.porosityMinimum,
            d.porosityAverage,
            d.porosityMaximum);

        printPhase_("Oil  ", Indices::Phase::liquid, d);
        printPhase_("Gas  ", Indices::Phase::vapor, d);
        if constexpr (Indices::hasWater)
            printPhase_("Water", Indices::Phase::water, d);

        if constexpr (Indices::fullyCompositionalThreePhase)
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "  Phase-state cells [O/G/W/OG/OW/GW/OGW]: %lld / %lld / %lld / %lld / %lld / %lld / %lld\n"
                "  Hysteresis-suppressed [O/G/W]       : %lld / %lld / %lld\n",
                d.phasePresenceCells[1], d.phasePresenceCells[2], d.phasePresenceCells[4],
                d.phasePresenceCells[3], d.phasePresenceCells[5], d.phasePresenceCells[6],
                d.phasePresenceCells[7],
                d.phaseSuppressionCells[0], d.phaseSuppressionCells[1],
                d.phaseSuppressionCells[2]);
        }
        else
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "  Hydrocarbon state cells [L/V/L+V]  : %lld / %lld / %lld\n",
                d.liquidOnlyCells,
                d.vaporOnlyCells,
                d.twoPhaseCells);
        }

        if constexpr (Indices::hasAqueousCO2Dissolution)
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "  Aqueous CO2 mass fraction [avg/max] : %.6e / %.6e\n",
                d.aqueousCO2MassFractionAverage,
                d.aqueousCO2MassFractionMaximum);
        }
        if constexpr (Indices::hasLandTrapping)
        {
            PetscPrintf(
                PETSC_COMM_SELF,
                "  Land trapped Sg [avg/max]            : %.6e / %.6e\n",
                d.trappedGasSaturationAverage,
                d.trappedGasSaturationMaximum);
        }

        PetscPrintf(PETSC_COMM_SELF,
                    "%s\n", consoleRule().c_str());
    }

    template <class Diagnostics>
    void printPhase_(const char *name, int phaseIndex, const Diagnostics &d) const
    {
        const std::size_t p = static_cast<std::size_t>(phaseIndex);
        PetscPrintf(
            PETSC_COMM_SELF,
            "  %s S[min/avg/max]=%.4g/%.4g/%.4g  rho=%.4g/%.4g/%.4g kg/m3  mu=%.3e/%.3e/%.3e Pa.s\n",
            name,
            d.saturationMinimum[p], d.saturationAverage[p], d.saturationMaximum[p],
            d.densityMinimum[p], d.densityAverage[p], d.densityMaximum[p],
            d.viscosityMinimum[p], d.viscosityAverage[p], d.viscosityMaximum[p]);
    }

    template <class Diagnostics>
    void append_(std::size_t step, double timeSeconds, const Diagnostics &d) const
    {
        if (rank_ != 0)
            return;

        const auto file = options_.resultDirectory / options_.historyFile;
        std::ofstream stream(file, std::ios::out | std::ios::app);
        if (!stream)
            throw std::runtime_error("Failed to append reservoir diagnostics file: " + file.string());
        stream << std::setprecision(16);

        const auto sat = [&](int phase, int which) {
            const auto p = static_cast<std::size_t>(phase);
            return which == 0 ? d.saturationMinimum[p] : (which == 1 ? d.saturationAverage[p] : d.saturationMaximum[p]);
        };
        const auto rho = [&](int phase, int which) {
            const auto p = static_cast<std::size_t>(phase);
            return which == 0 ? d.densityMinimum[p] : (which == 1 ? d.densityAverage[p] : d.densityMaximum[p]);
        };
        const auto mu = [&](int phase, int which) {
            const auto p = static_cast<std::size_t>(phase);
            return which == 0 ? d.viscosityMinimum[p] : (which == 1 ? d.viscosityAverage[p] : d.viscosityMaximum[p]);
        };
        double swMin = 0.0, swAvg = 0.0, swMax = 0.0;
        double rhoWMin = 0.0, rhoWAvg = 0.0, rhoWMax = 0.0;
        double muWMin = 0.0, muWAvg = 0.0, muWMax = 0.0;
        if constexpr (Indices::hasWater)
        {
            swMin = sat(Indices::Phase::water, 0);
            swAvg = sat(Indices::Phase::water, 1);
            swMax = sat(Indices::Phase::water, 2);
            rhoWMin = rho(Indices::Phase::water, 0);
            rhoWAvg = rho(Indices::Phase::water, 1);
            rhoWMax = rho(Indices::Phase::water, 2);
            muWMin = mu(Indices::Phase::water, 0);
            muWAvg = mu(Indices::Phase::water, 1);
            muWMax = mu(Indices::Phase::water, 2);
        }

        stream
            << step << ',' << timeSeconds / options_.secondsPerDisplayUnit << ','
            << d.cellCount << ',' << d.bulkVolume << ',' << d.poreVolume << ','
            << d.pressureMinimum / options_.pressureScale << ','
            << d.pressureAverage / options_.pressureScale << ','
            << d.pressureMaximum / options_.pressureScale << ','
            << d.porosityMinimum << ',' << d.porosityAverage << ',' << d.porosityMaximum << ','
            << sat(Indices::Phase::liquid,0) << ',' << sat(Indices::Phase::liquid,1) << ',' << sat(Indices::Phase::liquid,2) << ','
            << sat(Indices::Phase::vapor,0) << ',' << sat(Indices::Phase::vapor,1) << ',' << sat(Indices::Phase::vapor,2) << ','
            << swMin << ',' << swAvg << ',' << swMax << ','
            << rho(Indices::Phase::liquid,0) << ',' << rho(Indices::Phase::liquid,1) << ',' << rho(Indices::Phase::liquid,2) << ','
            << rho(Indices::Phase::vapor,0) << ',' << rho(Indices::Phase::vapor,1) << ',' << rho(Indices::Phase::vapor,2) << ','
            << rhoWMin << ',' << rhoWAvg << ',' << rhoWMax << ','
            << mu(Indices::Phase::liquid,0) << ',' << mu(Indices::Phase::liquid,1) << ',' << mu(Indices::Phase::liquid,2) << ','
            << mu(Indices::Phase::vapor,0) << ',' << mu(Indices::Phase::vapor,1) << ',' << mu(Indices::Phase::vapor,2) << ','
            << muWMin << ',' << muWAvg << ',' << muWMax << ',';

        if constexpr (Indices::fullyCompositionalThreePhase)
        {
            stream
                << d.phasePresenceCells[1] << ',' << d.phasePresenceCells[2] << ','
                << d.phasePresenceCells[4] << ',' << d.phasePresenceCells[3] << ','
                << d.phasePresenceCells[5] << ',' << d.phasePresenceCells[6] << ','
                << d.phasePresenceCells[7] << ','
                << d.phaseSuppressionCells[0] << ',' << d.phaseSuppressionCells[1] << ','
                << d.phaseSuppressionCells[2] << ',';
        }
        else
        {
            stream << d.liquidOnlyCells << ',' << d.vaporOnlyCells << ',' << d.twoPhaseCells << ',';
        }

        stream
            << d.aqueousCO2MassFractionAverage << ',' << d.aqueousCO2MassFractionMaximum << ','
            << d.trappedGasSaturationAverage << ',' << d.trappedGasSaturationMaximum << '\n';
    }

    ReservoirDiagnosticsOutputOptions options_;
    PetscMPIInt rank_{0};
};

} // namespace MPMC
