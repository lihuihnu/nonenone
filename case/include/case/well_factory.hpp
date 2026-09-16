/**
 * @file well_factory.hpp
 * @brief 统一正式算例的 Natural 井构建、控制映射和完井装配公共逻辑。
 *
 * The helper deliberately owns no case parameters.  Individual well_config.hpp
 * files remain the single source of physical inputs; this layer only translates
 * those definitions into NaturalWell objects and grid-specific perforations.
 */
#pragma once

#include <well/well.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace MPMC::cases
{

namespace detail
{

template <class T, class = void>
struct HasMinimumBhp : std::false_type {};
template <class T>
struct HasMinimumBhp<T, std::void_t<decltype(std::declval<const T &>().minimumBhp)>>
    : std::true_type {};

template <class T, class = void>
struct HasMaximumBhp : std::false_type {};
template <class T>
struct HasMaximumBhp<T, std::void_t<decltype(std::declval<const T &>().maximumBhp)>>
    : std::true_type {};

template <class T, class = void>
struct HasMaximumWaterRate : std::false_type {};
template <class T>
struct HasMaximumWaterRate<T, std::void_t<decltype(std::declval<const T &>().maximumWaterRate)>>
    : std::true_type {};

template <class T, class = void>
struct HasSchedule : std::false_type {};
template <class T>
struct HasSchedule<T,
                   std::void_t<decltype(std::declval<const T &>().scheduleEnabled),
                               decltype(std::declval<const T &>().openTime),
                               decltype(std::declval<const T &>().closeTime)>>
    : std::true_type {};

template <class T, class = void>
struct HasInjectionDefinition : std::false_type {};
template <class T>
struct HasInjectionDefinition<T,
                              std::void_t<decltype(std::declval<const T &>().injectionPhase),
                                          decltype(std::declval<const T &>().injectedComponent)>>
    : std::true_type {};

template <class Enum, class = void>
struct HasReservoirTotalRate : std::false_type {};
template <class Enum>
struct HasReservoirTotalRate<
    Enum, std::void_t<decltype(Enum::ReservoirTotalRate)>>
    : std::true_type {};

template <class Enum>
[[nodiscard]] WellControl toNaturalControl(Enum control)
{
    if constexpr (HasReservoirTotalRate<Enum>::value)
    {
        if (control == Enum::ReservoirTotalRate)
            return WellControl::ReservoirTotalRate;
    }
    switch (control)
    {
    case Enum::Bhp:       return WellControl::Bhp;
    case Enum::TotalRate: return WellControl::TotalRate;
    case Enum::OilRate:   return WellControl::OilRate;
    case Enum::GasRate:   return WellControl::GasRate;
    case Enum::WaterRate: return WellControl::WaterRate;
    default: break;
    }
    throw std::logic_error("Unknown case well control.");
}

template <class Enum>
[[nodiscard]] WellType toNaturalType(Enum type)
{
    switch (type)
    {
    case Enum::Injector: return WellType::Injector;
    case Enum::Producer: return WellType::Producer;
    }
    throw std::logic_error("Unknown case well type.");
}

template <class Indices, class Enum>
[[nodiscard]] int toNaturalPhase(Enum phase)
{
    switch (phase)
    {
    case Enum::Oil:   return Indices::Phase::liquid;
    case Enum::Gas:   return Indices::Phase::vapor;
    case Enum::Water: return Indices::Phase::water;
    }
    throw std::logic_error("Unknown case injection phase.");
}

template <class Indices, class CellId, class Definition>
void applyCommonDefinition(WellSpecification<Indices, CellId> &well, const Definition &definition)
{
    well.radius = definition.radius;
    well.skin = definition.skin;

    if constexpr (HasSchedule<Definition>::value)
    {
        well.schedule.enabled = definition.scheduleEnabled;
        well.schedule.openTime = definition.openTime;
        well.schedule.closeTime = definition.closeTime;
    }

    if constexpr (HasMinimumBhp<Definition>::value)
    {
        if (std::isfinite(definition.minimumBhp))
            well.setMinimumBhp(definition.minimumBhp);
    }
    if constexpr (HasMaximumBhp<Definition>::value)
    {
        if (std::isfinite(definition.maximumBhp))
            well.setMaximumBhp(definition.maximumBhp);
    }
    if constexpr (HasMaximumWaterRate<Definition>::value)
    {
        if (std::isfinite(definition.maximumWaterRate))
            well.setMaximumWaterRate(definition.maximumWaterRate);
    }

    if constexpr (HasInjectionDefinition<Definition>::value)
    {
        if (definition.injectedComponent >= 0)
        {
            const int phase = toNaturalPhase<Indices>(definition.injectionPhase);
            const int component = definition.injectedComponent;
            if (component >= Indices::numComponents)
                throw std::invalid_argument("Injected component index is outside the fluid system.");
            well.injectionPhaseFraction[static_cast<std::size_t>(phase)] = 1.0;
            well.injectionComponentMassFraction[static_cast<std::size_t>(component)] = 1.0;
        }
    }
}

} // namespace detail

/** Resolve -1 style structured-grid indices without duplicating case helpers. */
[[nodiscard]] inline int resolveStructuredIndex(int value, int count)
{
    const int resolved = value < 0 ? count + value : value;
    if (resolved < 0 || resolved >= count)
        throw std::out_of_range("Structured well completion index is outside the grid.");
    return resolved;
}

/**
 * Build one NaturalWell from a case definition after the caller has supplied
 * grid-specific perforations and a representative BHP cell.
 */
template <class Indices, class CellId, class Definition>
[[nodiscard]] WellSpecification<Indices, CellId> makeWell(
    const Definition &definition,
    CellId representativeCell,
    std::vector<WellPerforation<CellId>> perforations)
{
    WellSpecification<Indices, CellId> well(
        definition.id,
        definition.name,
        detail::toNaturalType(definition.type),
        detail::toNaturalControl(definition.control),
        definition.target,
        representativeCell,
        std::move(perforations));
    detail::applyCommonDefinition<Indices, CellId>(well, definition);
    return well;
}

/**
 * Build the common structured-grid vertical-well layout used by the reservoir
 * examples.  The callback supplies only the case-specific Peaceman WI because
 * rock properties may be constant or spatially varying.
 */
template <class Indices, class CellId, class Grid, class Definitions, class WellIndexFunction>
[[nodiscard]] std::vector<WellSpecification<Indices, CellId>> makeStructuredWells(
    Grid &grid,
    const Definitions &definitions,
    WellIndexFunction &&wellIndex)
{
    if constexpr (!Indices::hasWellUnknown)
        return {};

    const auto dims = grid.dimensions();
    std::vector<WellSpecification<Indices, CellId>> result;
    result.reserve(definitions.size());

    for (const auto &definition : definitions)
    {
        const int i = resolveStructuredIndex(definition.completion.i, dims[0]);
        const int j = resolveStructuredIndex(definition.completion.j, dims[1]);
        if (definition.completion.kBegin < 0 || definition.completion.kBegin >= dims[2])
            throw std::out_of_range("Structured well first completion layer is outside the grid.");
        const int count = definition.completion.kCount < 0
            ? dims[2] - definition.completion.kBegin
            : definition.completion.kCount;
        if (count <= 0 || definition.completion.kBegin + count > dims[2])
            throw std::out_of_range("Structured well completion interval is outside the grid.");

        const auto cellId = [&grid](int ii, int jj, int kk) {
            return static_cast<CellId>(grid.ijkToGlobal(ii, jj, kk));
        };
        auto perforations = makeVerticalPerforations<CellId>(
            i,
            j,
            definition.completion.kBegin,
            count,
            cellId,
            [&](int ii, int jj, int kk) {
                return wellIndex(grid, definition, ii, jj, kk);
            });

        result.push_back(makeWell<Indices, CellId>(
            definition,
            cellId(i, j, definition.completion.kBegin),
            std::move(perforations)));
    }
    return result;
}

} // namespace MPMC::cases
