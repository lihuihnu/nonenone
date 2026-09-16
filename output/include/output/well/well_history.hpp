/**
 * @file well_history.hpp
 * @brief 统一井历史 CSV 数据结构与写入接口。
 */
#pragma once

#include <output/core/writer.hpp>

#include <well/manager.hpp>
#include <well/state.hpp>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace MPMC
{

/**
 * @brief 将 Natural 井的各相流量写入 `well_phase_rates.csv`。
 *
 * states must be aligned with manager.wells().  The output order remains:
 *
 *     requested well id -> phase 0 -> phase 1 -> ...
 *
 * using the same phase-rate ordering as OutputWriter.
 */
template <
    class Indices,
    class Well>
void writeWellStates(
    OutputWriter<Indices> &writer,
    const WellManager<Well> &manager,
    const std::vector<WellState<Indices>> &states,
    const std::vector<int> &wellIds)
{
    if (states.size() != manager.size())
    {
        throw std::invalid_argument(
            "Well output states must be aligned with WellManager::wells().");
    }

    const auto &wells =
        manager.wells();

    for (int id : wellIds)
    {
        bool found = false;

        for (std::size_t index = 0;
             index < wells.size();
             ++index)
        {
            if (wells[index].id != id)
                continue;

            writer.writeWellPhaseRates(
                states[index].surfacePhaseRate);

            found = true;
            break;
        }

        if (!found)
        {
            throw std::out_of_range(
                "Requested well id does not exist in WellManager.");
        }
    }
}

} // namespace MPMC
