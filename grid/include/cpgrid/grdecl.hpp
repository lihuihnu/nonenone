/**
 * @file grdecl.hpp
 * @brief GRDECL 网格与岩石属性读取和标准化。
 */
#pragma once

#include <array>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace MPMC
{

/** Eclipse/GRDECL unit system relevant to grid geometry and permeability. */
enum class GrdeclUnitSystem
{
    Metric,
    Field,
    Lab
};

/** Options for direct GRDECL loading. */
struct GrdeclLoadOptions final
{
    /** Override a unit keyword found in the file. Empty means auto-detect. */
    std::optional<GrdeclUnitSystem> unitSystem;

    /** Convert coordinates to metre and permeability to m^2 for MPMC internals. */
    bool convertToSi{true};

    /** Geometry-node merging tolerance [m after SI conversion]. */
    double nodeMergeTolerance{1.0e-9};

    /** Keep only the largest Cartesian face-connected ACTNUM component. */
    bool keepLargestConnectedComponent{false};
};

struct GrdeclPoint final
{
    double x{0.0};
    double y{0.0};
    double z{0.0};
};

/** One active corner-point cell in canonical Cartesian ordering. */
struct GrdeclCell final
{
    int cartesianId{-1};
    int i{0};
    int j{0};
    int k{0};

    /** Corner index = ix + 2*iy + 4*iz. */
    std::array<GrdeclPoint, 8> corner{};

    /** Active-cell storage neighbour for I-, I+, J-, J+, K-, K+; -1 = boundary. */
    std::array<int, 6> neighborStorage{{-1, -1, -1, -1, -1, -1}};
};

/**
 * Parsed and geometry-expanded GRDECL input.
 *
 * Rock arrays are aligned with activeCells, not inactive Cartesian cells.
 * Geometry is converted to SI when GrdeclLoadOptions::convertToSi is true.
 */
struct GrdeclGridData final
{
    std::filesystem::path sourcePath;
    int nx{0};
    int ny{0};
    int nz{0};
    GrdeclUnitSystem unitSystem{GrdeclUnitSystem::Metric};
    bool unitSystemWasAssumed{false};
    bool convertedToSi{true};
    double nodeMergeTolerance{1.0e-9};

    std::size_t cartesianCellCount{0};
    std::size_t disconnectedCellCountRemoved{0};
    std::vector<int> actnum;
    std::vector<GrdeclCell> activeCells;

    std::vector<double> porosity;
    std::vector<double> permeabilityX;
    std::vector<double> permeabilityY;
    std::vector<double> permeabilityZ;

    [[nodiscard]] bool hasPorosity() const noexcept
    {
        return porosity.size() == activeCells.size();
    }

    [[nodiscard]] bool hasPermeability() const noexcept
    {
        return permeabilityX.size() == activeCells.size() &&
               permeabilityY.size() == activeCells.size() &&
               permeabilityZ.size() == activeCells.size();
    }
};

/** Parse an Eclipse-style GRDECL file, including simple INCLUDE files. */
[[nodiscard]] GrdeclGridData loadGrdecl(
    const std::filesystem::path &filename,
    const GrdeclLoadOptions &options = {});

[[nodiscard]] const char *grdeclUnitSystemName(GrdeclUnitSystem value) noexcept;

} // namespace MPMC
