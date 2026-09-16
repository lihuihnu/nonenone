/**
 * @file well_config.hpp
 * @brief 三维 SCW/CO2 共注井与上部生产井配置。
 */
#pragma once

#include "benchmark_common.hpp"

#include <array>

namespace WellConfig
{
using Type = ScwMigration3D::Type;
using Control = ScwMigration3D::Control;
using InjectionPhase = ScwMigration3D::InjectionPhase;
using Completion = ScwMigration3D::Completion;

struct WellDefinition : ScwMigration3D::BaseWellDefinition
{
    // The runtime fills the mixed SCW/CO2 stream explicitly.
    int injectedComponent{-1};
};

inline static constexpr std::array<WellDefinition, 2> wells{{
    {ScwMigration3D::baseWells[0], -1},
    {ScwMigration3D::baseWells[1], -1}
}};

} // namespace WellConfig
