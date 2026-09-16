/**
 * @file failure_diagnostics.hpp
 * @brief 失败 Newton 的近井组分质量守恒、面通量和井源项分解诊断。
 */
#pragma once

#include <adaptive_timestep/core/solve_result.hpp>
#include <common/units.hpp>
#include <well/types.hpp>

#include <petscsnes.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace MPMC::cases
{

/**
 * @brief rejected attempt 回滚前执行的质量守恒诊断 hook。
 *
 * 只诊断 MPI 全局绝对值最大的组分质量守恒行。旧的 phase-presence / flash
 * 诊断已经移除；当前版本直接把守恒式拆成 accumulation + face flux - well
 * source，并逐面逐相打印 potential/upwind/mobility/component flux。
 */
template <class Indices, class Runtime, class Config>
class FailedSolveDiagnostics final
{
public:
    FailedSolveDiagnostics(
        PetscMPIInt rank,
        std::string resultDirectory,
        bool print,
        bool write)
        : rank_(rank),
          resultDirectory_(std::move(resultDirectory)),
          print_(print),
          write_(write)
    {
        if (rank_ == 0 && write_)
        {
            std::filesystem::create_directories(resultDirectory_);
            initializeMassFile_();
            initializeFaceFile_();
            initializeWellFile_();
        }
    }

    void operator()(
        Runtime &runtime,
        SNES snes,
        Vec solution,
        double startTime,
        double timeStep,
        const MPMC::NonlinearSolveResult &result)
    {
        Vec residual = nullptr;
        PetscCallAbort(PETSC_COMM_WORLD, SNESGetFunction(snes, &residual, nullptr, nullptr));
        if (residual == nullptr)
            return;

        const auto diagnostic =
            runtime.evaluateResidualFailureDiagnostic(solution, residual);
        if (!diagnostic.valid)
            return;

        ++failureSerial_;
        if (rank_ == 0 && print_)
            printDiagnostic_(runtime, diagnostic, startTime, timeStep, result);
        if (rank_ == 0 && write_)
            writeDiagnostic_(runtime, diagnostic, startTime, timeStep, result);
    }

private:
    [[nodiscard]] static const char *phaseName_(int phase)
    {
        if (phase == Indices::Phase::liquid)
            return "O";
        if (phase == Indices::Phase::vapor)
            return "G";
        if constexpr (Indices::hasWater)
        {
            if (phase == Indices::Phase::water)
                return "W";
        }
        return "?";
    }

    [[nodiscard]] static const char *upwindName_(int side)
    {
        return side == 0 ? "SELF" : "NEIGHBOR";
    }

    [[nodiscard]] static std::string wellName_(const Runtime &runtime, int vectorIndex)
    {
        if (vectorIndex < 0 ||
            static_cast<std::size_t>(vectorIndex) >= runtime.wells().size())
            return "unknown";
        return runtime.wells()[static_cast<std::size_t>(vectorIndex)].name;
    }

    template <class Diagnostic>
    void printDiagnostic_(
        const Runtime &runtime,
        const Diagnostic &d,
        double startTime,
        double timeStep,
        const MPMC::NonlinearSolveResult &result) const
    {
        constexpr double day = 86400.0;
        constexpr double bar = MPMC::units::bar;

        PetscPrintf(PETSC_COMM_SELF,
                    "\n===================== FAILED MASS-BALANCE DIAGNOSTIC #%zu =====================\n",
                    failureSerial_);
        PetscPrintf(PETSC_COMM_SELF,
                    "  interval / dt                : %.12g -> %.12g day / %.12g day\n",
                    startTime / day,
                    (startTime + timeStep) / day,
                    timeStep / day);
        PetscPrintf(PETSC_COMM_SELF,
                    "  SNES reason                  : %s\n",
                    result.reason.empty() ? "unknown" : result.reason.c_str());
        PetscPrintf(PETSC_COMM_SELF,
                    "  max component |R|           : %.12e  component=%s\n",
                    d.maximumAbsoluteMassResidual,
                    Config::Fluid::componentNames[static_cast<std::size_t>(d.componentIndex)]);
        PetscPrintf(PETSC_COMM_SELF,
                    "  cell current/input           : %lld / %lld\n",
                    static_cast<long long>(d.currentCellId),
                    static_cast<long long>(d.inputCellId));
        PetscPrintf(PETSC_COMM_SELF,
                    "  volume / pressure            : %.12e m3 / %.9g bar\n",
                    d.cellVolume,
                    d.pressure / bar);

        PetscPrintf(PETSC_COMM_SELF,
                    "------------------------------------------------------------------------------\n");
        PetscPrintf(PETSC_COMM_SELF,
                    "  PHASE STATE USED BY FLUX\n");
        for (int phase = 0; phase < Indices::numPhases; ++phase)
        {
            const std::size_t p = static_cast<std::size_t>(phase);
            PetscPrintf(PETSC_COMM_SELF,
                        "    %s: S=%.9g rho=%.9e kg/m3 mobility=%.9e\n",
                        phaseName_(phase),
                        d.saturation[p],
                        d.density[p],
                        d.mobility[p]);
            std::ostringstream fractions;
            fractions << "       mass fractions:" << std::scientific << std::setprecision(6);
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                fractions << ' ' << Config::Fluid::componentNames[c]
                          << '=' << d.massFraction[p][c];
            }
            fractions << '\n';
            PetscPrintf(PETSC_COMM_SELF, "%s", fractions.str().c_str());
        }

        PetscPrintf(PETSC_COMM_SELF,
                    "------------------------------------------------------------------------------\n");
        PetscPrintf(PETSC_COMM_SELF,
                    "  COMPONENT BALANCE: R = ACCUMULATION + FACE_FLUX - WELL_SOURCE\n");
        PetscPrintf(PETSC_COMM_SELF,
                    "  %-10s %13s %13s %13s %13s %13s %13s\n",
                    "COMP", "ACC(kg/s)", "FACE(kg/s)", "WELL(kg/s)",
                    "RECON", "PETSC_R", "RECON-PETSC");
        for (int component = 0; component < Indices::numComponents; ++component)
        {
            const std::size_t c = static_cast<std::size_t>(component);
            PetscPrintf(PETSC_COMM_SELF,
                        "  %-10s %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
                        Config::Fluid::componentNames[c],
                        d.accumulationTerm[c],
                        d.faceFluxTerm[c],
                        d.wellSourceTerm[c],
                        d.reconstructedResidual[c],
                        d.petscResidual[c],
                        d.reconstructionError[c]);
            PetscPrintf(PETSC_COMM_SELF,
                        "       accumulation density current/previous = %.12e / %.12e kg/m3 bulk\n",
                        d.currentAccumulationDensity[c],
                        d.previousAccumulationDensity[c]);
        }

        PetscPrintf(PETSC_COMM_SELF,
                    "------------------------------------------------------------------------------\n");
        PetscPrintf(PETSC_COMM_SELF,
                    "  FACE-BY-FACE TPFA / UPWIND DETAILS  (positive = out of diagnostic cell)\n");
        for (std::size_t faceIndex = 0; faceIndex < d.faces.size(); ++faceIndex)
        {
            const auto &face = d.faces[faceIndex];
            PetscPrintf(PETSC_COMM_SELF,
                        "  FACE %zu: neighbor current/input=%lld/%lld  T=%.12e  gravityDz=%.12e\n",
                        faceIndex,
                        static_cast<long long>(face.neighborCellId),
                        static_cast<long long>(face.neighborInputCellId),
                        face.transmissibility,
                        face.gravityTerm);
            PetscPrintf(PETSC_COMM_SELF,
                        "          neighbor P=%.9g bar  S=",
                        face.neighborPressure / bar);
            for (int phase = 0; phase < Indices::numPhases; ++phase)
                PetscPrintf(PETSC_COMM_SELF, "%s%.9g%s",
                            phaseName_(phase),
                            face.neighborSaturation[static_cast<std::size_t>(phase)],
                            phase + 1 == Indices::numPhases ? "\n" : " ");

            for (int phase = 0; phase < Indices::numPhases; ++phase)
            {
                const std::size_t p = static_cast<std::size_t>(phase);
                PetscPrintf(PETSC_COMM_SELF,
                            "    phase %s: dPhi=% .9e bar  upwind=%-8s  lambda=%.9e  rho=%.9e  q=% .9e m3/s  mdot=% .9e kg/s\n",
                            phaseName_(phase),
                            face.potentialDifference[p] / bar,
                            upwindName_(face.upwindSide[p]),
                            face.upwindMobility[p],
                            face.upwindDensity[p],
                            face.darcyVolumeRate[p],
                            face.phaseMassRate[p]);
                std::ostringstream phaseFlux;
                phaseFlux << "       component flux:" << std::scientific << std::setprecision(6);
                for (int component = 0; component < Indices::numComponents; ++component)
                {
                    const std::size_t c = static_cast<std::size_t>(component);
                    phaseFlux << ' ' << Config::Fluid::componentNames[c]
                              << '=' << face.phaseComponentMassFlux[p][c];
                }
                phaseFlux << " kg/s\n";
                PetscPrintf(PETSC_COMM_SELF, "%s", phaseFlux.str().c_str());
            }

            std::ostringstream totalFlux;
            totalFlux << "    face total component flux:" << std::scientific << std::setprecision(6);
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                totalFlux << ' ' << Config::Fluid::componentNames[c]
                          << '=' << face.totalComponentMassFlux[c];
            }
            totalFlux << " kg/s\n";
            PetscPrintf(PETSC_COMM_SELF, "%s", totalFlux.str().c_str());
        }

        PetscPrintf(PETSC_COMM_SELF,
                    "------------------------------------------------------------------------------\n");
        PetscPrintf(PETSC_COMM_SELF,
                    "  WELL / PERFORATION SOURCE DETAILS  (positive = injection into reservoir)\n");
        if (d.perforations.empty())
        {
            PetscPrintf(PETSC_COMM_SELF, "    none on diagnostic cell\n");
        }
        for (std::size_t perfIndex = 0; perfIndex < d.perforations.size(); ++perfIndex)
        {
            const auto &perf = d.perforations[perfIndex];
            const auto type = static_cast<MPMC::WellType>(perf.wellType);
            const auto control = static_cast<MPMC::WellControl>(perf.wellControl);
            PetscPrintf(PETSC_COMM_SELF,
                        "  PERF %zu: well=%s id=%d type=%s control=%s WI=%.12e\n",
                        perfIndex,
                        wellName_(runtime, perf.wellVectorIndex).c_str(),
                        perf.wellId,
                        std::string(MPMC::wellTypeName(type)).c_str(),
                        std::string(MPMC::wellControlName(control)).c_str(),
                        perf.wellIndex);
            PetscPrintf(PETSC_COMM_SELF,
                        "          target=% .9e signed_target=% .9e actual=% .9e control_R=% .9e BHP=%.9g bar\n",
                        perf.target,
                        perf.signedTarget,
                        perf.controlledRate,
                        perf.controlResidual,
                        perf.bhp / bar);

            for (int phase = 0; phase < Indices::numPhases; ++phase)
            {
                const std::size_t p = static_cast<std::size_t>(phase);
                PetscPrintf(PETSC_COMM_SELF,
                            "    phase %s: injFrac=%.9g dP=% .9e bar lambda=%.9e rho=%.9e q_sc=% .9e q_res=% .9e mdot=% .9e\n",
                            phaseName_(phase),
                            perf.injectionPhaseFraction[p],
                            perf.pressureDrop[p] / bar,
                            perf.mobility[p],
                            perf.density[p],
                            perf.surfacePhaseRate[p],
                            perf.reservoirPhaseRate[p],
                            perf.phaseMassRate[p]);
            }

            std::ostringstream injComp;
            injComp << "       injection component mass fraction:" << std::scientific << std::setprecision(6);
            std::ostringstream sourceComp;
            sourceComp << "       perforation component source   :" << std::scientific << std::setprecision(6);
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                injComp << ' ' << Config::Fluid::componentNames[c]
                        << '=' << perf.injectionComponentMassFraction[c];
                sourceComp << ' ' << Config::Fluid::componentNames[c]
                           << '=' << perf.componentMassSource[c];
            }
            injComp << '\n';
            sourceComp << " kg/s\n";
            PetscPrintf(PETSC_COMM_SELF, "%s", injComp.str().c_str());
            PetscPrintf(PETSC_COMM_SELF, "%s", sourceComp.str().c_str());
        }

        PetscPrintf(PETSC_COMM_SELF,
                    "==============================================================================\n");
    }

    void initializeMassFile_() const
    {
        const auto path = resultDirectory_ / "failed_mass_balance_diagnostics.csv";
        std::ofstream stream(path, std::ios::out | std::ios::trunc);
        if (!stream)
            throw std::runtime_error("Failed to create mass-balance diagnostic file: " + path.string());
        stream << "failure,start_day,end_day,dt_day,snes_reason,cell_id,input_cell_id,target_component,component,"
                  "current_accum_density,previous_accum_density,accumulation_kg_s,face_flux_kg_s,"
                  "well_source_kg_s,reconstructed_residual,petsc_residual,reconstruction_error\n";
    }

    void initializeFaceFile_() const
    {
        const auto path = resultDirectory_ / "failed_face_flux_diagnostics.csv";
        std::ofstream stream(path, std::ios::out | std::ios::trunc);
        if (!stream)
            throw std::runtime_error("Failed to create face-flux diagnostic file: " + path.string());
        stream << "failure,cell_id,input_cell_id,face_index,neighbor_cell_id,neighbor_input_cell_id,"
                  "transmissibility,gravity_dz,phase,dphi_bar,upwind,mobility,density_kg_m3,"
                  "darcy_rate_m3_s,phase_mass_rate_kg_s,component,phase_component_flux_kg_s,"
                  "face_total_component_flux_kg_s\n";
    }

    void initializeWellFile_() const
    {
        const auto path = resultDirectory_ / "failed_well_source_diagnostics.csv";
        std::ofstream stream(path, std::ios::out | std::ios::trunc);
        if (!stream)
            throw std::runtime_error("Failed to create well-source diagnostic file: " + path.string());
        stream << "failure,cell_id,input_cell_id,perf_index,well_name,well_id,type,control,target,signed_target,"
                  "actual_control_value,control_residual,bhp_bar,WI,phase,injection_phase_fraction,dp_bar,"
                  "mobility,density_kg_m3,surface_rate_m3_s,reservoir_rate_m3_s,phase_mass_rate_kg_s,"
                  "component,injection_component_mass_fraction,component_source_kg_s\n";
    }

    template <class Diagnostic>
    void writeDiagnostic_(
        const Runtime &runtime,
        const Diagnostic &d,
        double startTime,
        double timeStep,
        const MPMC::NonlinearSolveResult &result) const
    {
        constexpr double day = 86400.0;
        constexpr double bar = MPMC::units::bar;

        {
            const auto path = resultDirectory_ / "failed_mass_balance_diagnostics.csv";
            std::ofstream stream(path, std::ios::out | std::ios::app);
            if (!stream)
                throw std::runtime_error("Failed to append mass-balance diagnostic file: " + path.string());
            stream << std::setprecision(16);
            for (int component = 0; component < Indices::numComponents; ++component)
            {
                const std::size_t c = static_cast<std::size_t>(component);
                stream << failureSerial_ << ','
                       << startTime / day << ','
                       << (startTime + timeStep) / day << ','
                       << timeStep / day << ','
                       << result.reason << ','
                       << d.currentCellId << ','
                       << d.inputCellId << ','
                       << Config::Fluid::componentNames[static_cast<std::size_t>(d.componentIndex)] << ','
                       << Config::Fluid::componentNames[c] << ','
                       << d.currentAccumulationDensity[c] << ','
                       << d.previousAccumulationDensity[c] << ','
                       << d.accumulationTerm[c] << ','
                       << d.faceFluxTerm[c] << ','
                       << d.wellSourceTerm[c] << ','
                       << d.reconstructedResidual[c] << ','
                       << d.petscResidual[c] << ','
                       << d.reconstructionError[c] << '\n';
            }
        }

        {
            const auto path = resultDirectory_ / "failed_face_flux_diagnostics.csv";
            std::ofstream stream(path, std::ios::out | std::ios::app);
            if (!stream)
                throw std::runtime_error("Failed to append face-flux diagnostic file: " + path.string());
            stream << std::setprecision(16);
            for (std::size_t faceIndex = 0; faceIndex < d.faces.size(); ++faceIndex)
            {
                const auto &face = d.faces[faceIndex];
                for (int phase = 0; phase < Indices::numPhases; ++phase)
                {
                    const std::size_t p = static_cast<std::size_t>(phase);
                    for (int component = 0; component < Indices::numComponents; ++component)
                    {
                        const std::size_t c = static_cast<std::size_t>(component);
                        stream << failureSerial_ << ','
                               << d.currentCellId << ',' << d.inputCellId << ','
                               << faceIndex << ',' << face.neighborCellId << ','
                               << face.neighborInputCellId << ','
                               << face.transmissibility << ',' << face.gravityTerm << ','
                               << phaseName_(phase) << ',' << face.potentialDifference[p] / bar << ','
                               << upwindName_(face.upwindSide[p]) << ','
                               << face.upwindMobility[p] << ',' << face.upwindDensity[p] << ','
                               << face.darcyVolumeRate[p] << ',' << face.phaseMassRate[p] << ','
                               << Config::Fluid::componentNames[c] << ','
                               << face.phaseComponentMassFlux[p][c] << ','
                               << face.totalComponentMassFlux[c] << '\n';
                    }
                }
            }
        }

        {
            const auto path = resultDirectory_ / "failed_well_source_diagnostics.csv";
            std::ofstream stream(path, std::ios::out | std::ios::app);
            if (!stream)
                throw std::runtime_error("Failed to append well-source diagnostic file: " + path.string());
            stream << std::setprecision(16);
            for (std::size_t perfIndex = 0; perfIndex < d.perforations.size(); ++perfIndex)
            {
                const auto &perf = d.perforations[perfIndex];
                const auto type = static_cast<MPMC::WellType>(perf.wellType);
                const auto control = static_cast<MPMC::WellControl>(perf.wellControl);
                for (int phase = 0; phase < Indices::numPhases; ++phase)
                {
                    const std::size_t p = static_cast<std::size_t>(phase);
                    for (int component = 0; component < Indices::numComponents; ++component)
                    {
                        const std::size_t c = static_cast<std::size_t>(component);
                        stream << failureSerial_ << ','
                               << d.currentCellId << ',' << d.inputCellId << ',' << perfIndex << ','
                               << wellName_(runtime, perf.wellVectorIndex) << ',' << perf.wellId << ','
                               << MPMC::wellTypeName(type) << ',' << MPMC::wellControlName(control) << ','
                               << perf.target << ',' << perf.signedTarget << ',' << perf.controlledRate << ','
                               << perf.controlResidual << ',' << perf.bhp / bar << ',' << perf.wellIndex << ','
                               << phaseName_(phase) << ',' << perf.injectionPhaseFraction[p] << ','
                               << perf.pressureDrop[p] / bar << ',' << perf.mobility[p] << ','
                               << perf.density[p] << ',' << perf.surfacePhaseRate[p] << ','
                               << perf.reservoirPhaseRate[p] << ',' << perf.phaseMassRate[p] << ','
                               << Config::Fluid::componentNames[c] << ','
                               << perf.injectionComponentMassFraction[c] << ','
                               << perf.componentMassSource[c] << '\n';
                    }
                }
            }
        }
    }

    PetscMPIInt rank_{0};
    std::filesystem::path resultDirectory_;
    bool print_{true};
    bool write_{true};
    std::size_t failureSerial_{0};
};

} // namespace MPMC::cases
