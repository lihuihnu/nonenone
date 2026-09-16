/**
 * @file cpgrid_saver.hpp
 * @brief CpGrid 结果按原始输入单元顺序保存的 PETSc 输出器。
 */
#pragma once

#include <output/core/writer.hpp>
#include <output/core/types.hpp>
#include <output/petsc/cpgrid_output_access.hpp>
#include <output/well/well_history.hpp>

#include <well/manager.hpp>
#include <well/state.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <filesystem>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief CpGrid/PETSc 的轻量 CSV 输出门面。
 *
 * 仅负责两件事：通过 CpGridOutputAccess 做 collective 采样/归约；通过
 * OutputWriter 在 rank 0 格式化写出。它不拥有 Grid、Vec、WellManager
 * 或模拟状态。
 */
template <class Indices, class Grid>
class CpGridSaver final
{
public:
    explicit CpGridSaver(
        std::filesystem::path outputDirectory,
        const Grid &grid)
        : access_(grid),
          writer_(
              std::move(outputDirectory),
              grid.mesh().rank() == 0)
    {
    }

    CpGridSaver(
        const CpGridSaver &) = delete;
    CpGridSaver &operator=(
        const CpGridSaver &) = delete;
    CpGridSaver(
        CpGridSaver &&) = delete;
    CpGridSaver &operator=(
        CpGridSaver &&) = delete;

    ~CpGridSaver() = default;

    void saveSolution(
        Vec solution,
        const std::vector<OutputCellSelection> &cells,
        const std::vector<PetscInt> &variables)
    {
        writer_.writeSolutionSamples(
            access_.sample(
                solution,
                static_cast<PetscInt>(
                    Indices::numPrimaryVariables),
                cells,
                variables));
    }

    void savePhaseState(
        Vec phaseState,
        const std::vector<OutputCellSelection> &cells,
        const std::vector<PetscInt> &variables)
    {
        writer_.writePhaseStateSamples(
            access_.sample(
                phaseState,
                static_cast<PetscInt>(
                    Indices::numPhaseStateVariables),
                cells,
                variables));
    }

    void saveMass(
        MassSeries series,
        Vec field,
        const std::vector<OutputCellSelection> &cells,
        const std::vector<PetscInt> &components)
    {
        requireMassFeature_(series);

        writer_.writeMassSamples(
            series,
            access_.sample(
                field,
                static_cast<PetscInt>(
                    Indices::numComponents),
                cells,
                components));
    }

    void saveMassTotal(
        MassSeries series,
        Vec field,
        double currentTime)
    {
        requireMassFeature_(series);

        writer_.writeMassTotals(
            series,
            currentTime,
            access_.template componentTotals<
                static_cast<std::size_t>(
                    Indices::numComponents)>(
                        field));
    }

    template <class Well>
    void saveWell(
        const WellManager<Well> &manager,
        const std::vector<
            WellState<Indices>> &states,
        const std::vector<int> &wellIds)
    {
        writeWellStates(
            writer_,
            manager,
            states,
            wellIds);
    }


private:
    static void requireMassFeature_(MassSeries series)
    {
        if (!massSeriesEnabled<Indices>(series))
            throw std::logic_error(
                "Requested mass output is disabled by the model configuration.");
    }

    CpGridOutputAccess<Grid> access_;
    OutputWriter<Indices> writer_;
};

} // namespace MPMC
