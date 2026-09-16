/**
 * @file mesh_grdecl_input.cpp
 * @brief GRDECL corner-point 数据到 CanonicalMeshData 的输入适配器。
 */
#include "mesh_canonical_data.hpp"

#include <cpgrid/grdecl.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace MPMC::detail
{

namespace
{

struct QuantizedPointKey final
{
    long long x{0};
    long long y{0};
    long long z{0};

    bool operator==(const QuantizedPointKey &other) const noexcept
    {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct QuantizedPointHash final
{
    std::size_t operator()(const QuantizedPointKey &key) const noexcept
    {
        std::size_t h = std::hash<long long>{}(key.x);
        h ^= std::hash<long long>{}(key.y) + 0x9e3779b97f4a7c15ULL + (h << 6U) + (h >> 2U);
        h ^= std::hash<long long>{}(key.z) + 0x9e3779b97f4a7c15ULL + (h << 6U) + (h >> 2U);
        return h;
    }
};

QuantizedPointKey quantize(
    const GrdeclPoint &point,
    double tolerance)
{
    const auto bucket = [tolerance](double coordinate) {
        const double scaled = std::floor(coordinate / tolerance);
        if (!std::isfinite(scaled) ||
            scaled <= static_cast<double>(std::numeric_limits<long long>::min()) ||
            scaled >= static_cast<double>(std::numeric_limits<long long>::max()))
        {
            throw std::overflow_error(
                "GRDECL coordinate exceeds the node-merge index range.");
        }
        return static_cast<long long>(scaled);
    };

    return QuantizedPointKey{
        bucket(point.x),
        bucket(point.y),
        bucket(point.z)};
}

bool isSameNode(
    const GrdeclPoint &point,
    const CanonicalMeshNode &node,
    double tolerance) noexcept
{
    const double dx = point.x - node.x;
    const double dy = point.y - node.y;
    const double dz = point.z - node.z;
    return dx * dx + dy * dy + dz * dz <= tolerance * tolerance;
}

} // namespace

CanonicalMeshData makeGrdeclCanonicalMeshData(
    const GrdeclGridData &grdecl)
{
    if (grdecl.activeCells.empty())
        throw std::runtime_error(
            "GRDECL mesh contains no active cells.");
    if (!(grdecl.nodeMergeTolerance > 0.0) ||
        !std::isfinite(grdecl.nodeMergeTolerance))
    {
        throw std::invalid_argument(
            "GRDECL node merge tolerance is invalid.");
    }

    CanonicalMeshData result;
    result.sourcePath = grdecl.sourcePath.string();
    result.sourceFormat = "GRDECL";
    result.logicalDimensions = {{grdecl.nx, grdecl.ny, grdecl.nz}};

    const std::size_t cellCount = grdecl.activeCells.size();
    result.cells.reserve(cellCount);

    std::unordered_multimap<
        QuantizedPointKey,
        std::size_t,
        QuantizedPointHash>
        nodeByKey;
    nodeByKey.reserve(cellCount * 8);

    std::vector<std::array<std::size_t, 8>> nodeIndexByCell(cellCount);
    for (std::size_t storage = 0; storage < cellCount; ++storage)
    {
        const auto &input = grdecl.activeCells[storage];
        for (std::size_t c = 0; c < input.corner.size(); ++c)
        {
            const auto key =
                quantize(input.corner[c], grdecl.nodeMergeTolerance);
            std::size_t nodeIndex =
                std::numeric_limits<std::size_t>::max();

            for (long long dz = -1; dz <= 1; ++dz)
                for (long long dy = -1; dy <= 1; ++dy)
                    for (long long dx = -1; dx <= 1; ++dx)
                    {
                        const QuantizedPointKey candidate{
                            key.x + dx,
                            key.y + dy,
                            key.z + dz};
                        const auto range = nodeByKey.equal_range(candidate);
                        for (auto found = range.first;
                             found != range.second;
                             ++found)
                        {
                            if (isSameNode(
                                    input.corner[c],
                                    result.nodes.at(found->second),
                                    grdecl.nodeMergeTolerance))
                            {
                                nodeIndex = std::min(
                                    nodeIndex,
                                    found->second);
                            }
                        }
                    }

            if (nodeIndex == std::numeric_limits<std::size_t>::max())
            {
                const std::size_t index = result.nodes.size();
                const auto &point = input.corner[c];
                result.nodes.push_back(
                    CanonicalMeshNode{point.x, point.y, point.z});
                nodeByKey.emplace(key, index);
                nodeIndex = index;
            }

            nodeIndexByCell[storage][c] = nodeIndex;
        }
    }

    std::vector<std::array<CanonicalMeshIndex, 6>> faceId(cellCount);
    for (auto &ids : faceId)
        ids.fill(-1);

    constexpr std::array<int, 6> opposite{{1, 0, 3, 2, 5, 4}};
    for (std::size_t storage = 0; storage < cellCount; ++storage)
    {
        const auto &input = grdecl.activeCells[storage];
        for (std::size_t direction = 0; direction < 6; ++direction)
        {
            if (faceId[storage][direction] >= 0)
                continue;

            const int neighbor = input.neighborStorage[direction];
            const CanonicalMeshIndex id =
                static_cast<CanonicalMeshIndex>(result.faceNeighbors.size());
            faceId[storage][direction] = id;
            result.faceNeighbors.push_back(
                {{static_cast<CanonicalMeshIndex>(storage),
                  static_cast<CanonicalMeshIndex>(neighbor)}});

            if (neighbor >= 0)
            {
                const auto n = static_cast<std::size_t>(neighbor);
                if (n >= cellCount)
                {
                    throw std::runtime_error(
                        "GRDECL topology references an invalid active-cell storage index.");
                }
                faceId[n][static_cast<std::size_t>(opposite[direction])] = id;
            }
        }
    }

    constexpr std::array<std::array<int, 4>, 6> cornersByFace{{
        {{0, 2, 6, 4}},
        {{1, 5, 7, 3}},
        {{0, 4, 5, 1}},
        {{2, 3, 7, 6}},
        {{0, 1, 3, 2}},
        {{4, 6, 7, 5}}
    }};

    std::unordered_set<CanonicalMeshIndex> cartesianIds;
    cartesianIds.reserve(cellCount);
    for (std::size_t storage = 0; storage < cellCount; ++storage)
    {
        const auto &input = grdecl.activeCells[storage];
        CanonicalMeshCell cell;
        cell.cartesianId = static_cast<CanonicalMeshIndex>(input.cartesianId);
        if (!cartesianIds.insert(cell.cartesianId).second)
            throw std::runtime_error(
                "GRDECL contains duplicate Cartesian cell id.");
        cell.faces.reserve(6);

        for (std::size_t direction = 0; direction < 6; ++direction)
        {
            CanonicalMeshFace face;
            face.inputPosition = static_cast<CanonicalMeshIndex>(direction + 1);
            face.inputFaceId = faceId[storage][direction];
            face.nodeIndices.reserve(4);
            for (int localCorner : cornersByFace[direction])
            {
                face.nodeIndices.push_back(
                    nodeIndexByCell[storage][static_cast<std::size_t>(localCorner)]);
            }
            cell.faces.push_back(std::move(face));
        }

        result.cells.push_back(std::move(cell));
    }

    return result;
}

} // namespace MPMC::detail
