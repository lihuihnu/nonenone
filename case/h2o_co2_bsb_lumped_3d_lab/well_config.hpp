/**
 * @file well_config.hpp
 * @brief 定义实验室砂箱算例的缩放布井、完井区间和井控条件。
 */
#pragma once

#include "case_config.hpp"

#include <case/well_config_schema.hpp>

#include <array>
#include <stdexcept>

namespace LabScw3DWell
{

using Type = MPMC::cases::well_config::Type;
using Control = MPMC::cases::well_config::Control;
using InjectionPhase = MPMC::cases::well_config::InjectionPhase;
using WellDefinition =
    MPMC::cases::well_config::ScheduledStructuredWellDefinition<LabScw3D::Config>;

inline std::array<WellDefinition, 2> definitions(
    int nx, int ny, int nz, double totalReservoirRate)
{
    if (nx != LabScw3D::Grid::nx || ny != LabScw3D::Grid::ny ||
        nz != LabScw3D::Grid::nz)
        throw std::invalid_argument(
            "This comparison must use the unchanged H2O-CO2-nC10 60x20x1 grid.");

    std::array<WellDefinition, 2> result{};
    result[0] = {0, "MIXED_INJ", Type::Injector, Control::TotalRate,
                 totalReservoirRate, LabScw3D::InitialState::pressure,
                 {0, 9, 0, 1}, 0.10, 0.0,
                 InjectionPhase::Water, -1};
    result[0].maximumBhp = BenchmarkCommon::injectorMaximumBhp;

    result[1] = {1, "PROD", Type::Producer, Control::Bhp,
                 BenchmarkCommon::producerBhp, BenchmarkCommon::producerBhp,
                 {59, 9, 0, 1}, 0.10, 0.0,
                 InjectionPhase::Oil, -1};
    return result;
}

} // namespace LabScw3DWell
