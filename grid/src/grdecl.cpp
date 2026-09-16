/**
 * @file grdecl.cpp
 * @brief GRDECL 公共加载流程的高层编排实现。
 */
#include <cpgrid/grdecl.hpp>

#include "grdecl_detail.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MPMC
{

const char *grdeclUnitSystemName(GrdeclUnitSystem value) noexcept
{
    switch (value)
    {
    case GrdeclUnitSystem::Metric: return "METRIC";
    case GrdeclUnitSystem::Field: return "FIELD";
    case GrdeclUnitSystem::Lab: return "LAB";
    }
    return "UNKNOWN";
}

GrdeclGridData loadGrdecl(
    const std::filesystem::path &filename,
    const GrdeclLoadOptions &options)
{
    using namespace grdecl_detail;

    if (!(options.nodeMergeTolerance > 0.0) || !std::isfinite(options.nodeMergeTolerance))
        throw std::invalid_argument("GRDECL nodeMergeTolerance must be finite and positive.");

    const RawDeck deck = parseDeck(filename);

    const auto dims = dimensions(deck);
    const int nx = dims[0];
    const int ny = dims[1];
    const int nz = dims[2];
    const std::size_t cartesianCount =
        static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz);
    if (cartesianCount > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error(
            "GRDECL Cartesian cell count exceeds the supported integer index range.");

    const std::size_t expectedCoord =
        6 * (static_cast<std::size_t>(nx) + 1) *
        (static_cast<std::size_t>(ny) + 1);
    const std::size_t expectedZcorn = 8 * cartesianCount;
    auto coord = expandReal(required(deck, "COORD"), "COORD", expectedCoord);
    auto zcorn = expandReal(required(deck, "ZCORN"), "ZCORN", expectedZcorn);

    GrdeclGridData result;
    result.sourcePath = std::filesystem::absolute(filename).lexically_normal();
    result.nx = nx;
    result.ny = ny;
    result.nz = nz;
    result.cartesianCellCount = cartesianCount;
    result.convertedToSi = options.convertToSi;
    result.nodeMergeTolerance = options.nodeMergeTolerance;

    if (options.unitSystem)
        result.unitSystem = *options.unitSystem;
    else if (const auto gridUnits = keywordLengthUnit(deck, "GRIDUNIT"))
        result.unitSystem = *gridUnits;
    else if (deck.units)
        result.unitSystem = *deck.units;
    else
    {
        result.unitSystem = GrdeclUnitSystem::Metric;
        result.unitSystemWasAssumed = true;
    }

    if (options.convertToSi)
        convertGeometryToSi(coord, zcorn, result.unitSystem);

    if (const auto mapAxes = optionalReal(deck, "MAPAXES", 6))
    {
        const auto mapUnits = keywordLengthUnit(deck, "MAPUNITS");
        const double mapFactor = mapUnits ? geometryLengthFactor(*mapUnits) : 1.0;
        const double gridFactor = geometryLengthFactor(result.unitSystem);
        const double originFactor = options.convertToSi ? mapFactor : mapFactor / gridFactor;
        applyMapAxes(coord, *mapAxes, originFactor);
    }

    const auto actnumIt = deck.keywordValues.find("ACTNUM");
    if (actnumIt == deck.keywordValues.end())
        result.actnum.assign(cartesianCount, 1);
    else
        result.actnum = expandInt(actnumIt->second, "ACTNUM", cartesianCount);

    if (options.keepLargestConnectedComponent)
        result.disconnectedCellCountRemoved =
            retainLargestConnectedComponent(result.actnum, nx, ny, nz);

    auto activeGrid = buildActiveGrid(coord, zcorn, result.actnum, nx, ny, nz);
    result.activeCells = std::move(activeGrid.activeCells);
    const auto &storageByCartesian = activeGrid.storageByCartesian;

    result.porosity = activeProperty(
        optionalReal(deck, "PORO", cartesianCount), storageByCartesian, cartesianCount, "PORO");
    result.permeabilityX = activeProperty(
        optionalReal(deck, "PERMX", cartesianCount), storageByCartesian, cartesianCount, "PERMX");
    result.permeabilityY = activeProperty(
        optionalReal(deck, "PERMY", cartesianCount), storageByCartesian, cartesianCount, "PERMY");
    result.permeabilityZ = activeProperty(
        optionalReal(deck, "PERMZ", cartesianCount), storageByCartesian, cartesianCount, "PERMZ");

    // 输入约定：缺失方向渗透率时沿用 PERMX，保持各向同性 deck 的旧行为。
    if (!result.permeabilityX.empty())
    {
        if (result.permeabilityY.empty())
            result.permeabilityY = result.permeabilityX;
        if (result.permeabilityZ.empty())
            result.permeabilityZ = result.permeabilityX;
    }

    if (options.convertToSi)
    {
        convertPermeabilityToSi(result.permeabilityX);
        convertPermeabilityToSi(result.permeabilityY);
        convertPermeabilityToSi(result.permeabilityZ);
    }

    for (double phi : result.porosity)
        if (!std::isfinite(phi) || phi < 0.0 || phi > 1.0)
            throw std::runtime_error("PORO must be finite and inside [0,1].");
    for (const auto *perm : {&result.permeabilityX, &result.permeabilityY, &result.permeabilityZ})
        for (double value : *perm)
            if (!std::isfinite(value) || value < 0.0)
                throw std::runtime_error("Permeability must be finite and non-negative.");

    return result;
}

} // namespace MPMC
