/**
 * @file types.hpp
 * @brief CpGrid 公共类型、索引和轻量别名。
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace MPMC
{

/** @brief 选取输出单元时使用的坐标/编号空间。 */
enum class OutputCellIdSpace
{
    /** Original zero-based row in the grid data files; canonical external id. */
    InputIndex,
    /** Partition/current id used internally by PETSc. */
    CurrentId,
    /** MRST Cartesian id from G.cells.indexMap. */
    CartesianId,
    /** Backward-compatible name; identical to InputIndex. */
    StorageIndex = InputIndex
};

/**
 * @brief 一个请求输出的单元及其编号约定。
 *
 * External sampling defaults to InputIndex so case files and post-processing
 * use the same numbering as the original grid/rock CSV rows.
 */
struct OutputCellSelection final
{
    std::int64_t id{0};
    OutputCellIdSpace space{OutputCellIdSpace::InputIndex};
};

/** @brief 模拟器可选写出的逐组分质量库存类型。 */
enum class MassSeries : std::size_t
{
    Adsorbed,
    Trapped,
    Dissolved,
    OilPhaseComponents,
    Count
};

/** @brief OutputWriter 与 PETSc 适配器共用的单个质量序列文件布局。 */
struct MassSeriesMetadata final
{
    std::string_view fieldFilename;
    std::string_view totalFilename;
    std::string_view totalColumn;
};

/**
 * @brief 所有质量库存序列的标准文件名与格式定义。
 *
 * Keeping this table in one place avoids filename/switch duplication between
 * rank-local writers, CpGrid savers and global inventory diagnostics.
 */
inline constexpr std::array<MassSeriesMetadata,
                            static_cast<std::size_t>(MassSeries::Count)>
    massSeriesMetadata{{
        {"adsorbed_mass_samples.csv", "adsorbed_mass_totals.csv", "total_adsorbed_mass_kg"},
        {"trapped_mass_samples.csv", "trapped_mass_totals.csv", "total_trapped_mass_kg"},
        {"dissolved_mass_samples.csv", "dissolved_mass_totals.csv", "total_dissolved_mass_kg"},
        {"oil_phase_component_mass_samples.csv", "oil_phase_component_mass_totals.csv",
         "total_oil_phase_component_mass_kg"},
    }};

/** @brief 将强类型质量序列 id 转换为数组索引。 */
[[nodiscard]] constexpr std::size_t massSeriesIndex(MassSeries series) noexcept
{
    return static_cast<std::size_t>(series);
}

/** @brief 返回指定质量库存序列的标准文件元数据。 */
[[nodiscard]] constexpr const MassSeriesMetadata &
massSeriesInfo(MassSeries series) noexcept
{
    return massSeriesMetadata[massSeriesIndex(series)];
}

/**
 * @brief 判断当前编译期模型特性是否存在指定质量库存。
 *
 * Oil-phase component mass is always meaningful. Adsorbed, trapped and
 * dissolved inventories are emitted only when the corresponding physical
 * model is enabled in `Indices`.
 */
template <class Indices>
[[nodiscard]] constexpr bool massSeriesEnabled(MassSeries series) noexcept
{
    switch (series)
    {
    case MassSeries::Adsorbed:
        return Indices::hasAdsorption;
    case MassSeries::Trapped:
        return Indices::hasLandTrapping;
    case MassSeries::Dissolved:
        return Indices::hasAqueousCO2Dissolution;
    case MassSeries::OilPhaseComponents:
        return true;
    case MassSeries::Count:
        return false;
    }
    return false;
}

} // namespace MPMC
