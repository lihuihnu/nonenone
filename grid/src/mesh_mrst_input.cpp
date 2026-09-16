/**
 * @file mesh_mrst_input.cpp
 * @brief MRST CSV 到 CanonicalMeshData 的输入适配器。
 */
#include "mesh_canonical_data.hpp"

#include <cpgrid/config.hpp>
#include <cpgrid/csv_reader.hpp>

#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace MPMC::detail
{

namespace
{

std::string path(
    const std::string &directory,
    const std::string &relative)
{
    return directory + "/" + relative;
}

} // namespace

CanonicalMeshData loadMrstCanonicalMeshData(
    const std::string &dataDirectory)
{
    CanonicalMeshData result;
    result.sourcePath = dataDirectory;
    result.sourceFormat = "MRST_CSV";

    const auto coordinates =
        readCsv<double>(
            path(dataDirectory, "G/nodes/coords/data.csv"));

    result.nodes.reserve(coordinates.size());
    for (const auto &row : coordinates)
    {
#if MESH_DIM == 3
        if (row.size() < 3)
            throw std::runtime_error(
                "3D node row has fewer than three coordinates.");

        result.nodes.push_back(
            CanonicalMeshNode{row[0], row[1], row[2]});
#else
        if (row.size() < 2)
            throw std::runtime_error(
                "2D node row has fewer than two coordinates.");

        result.nodes.push_back(
            CanonicalMeshNode{row[0], row[1], 0.0});
#endif
    }

    const auto indexMap =
        readMatlabIndexCsv<int>(
            path(dataDirectory, "G/cells/indexMap/data.csv"));
    const auto cellFaces =
        readMatlabIndexCsv<int>(
            path(dataDirectory, "G/cells/faces/data.csv"));
    const auto cellFacePos =
        readMatlabIndexCsv<int>(
            path(dataDirectory, "G/cells/facePos/data.csv"));
    const auto faceNodes =
        readMatlabIndexCsv<int>(
            path(dataDirectory, "G/faces/nodes/data.csv"));
    const auto faceNodePos =
        readMatlabIndexCsv<int>(
            path(dataDirectory, "G/faces/nodePos/data.csv"));

    if (cellFacePos.size() < indexMap.size() + 1)
        throw std::runtime_error(
            "G.cells.facePos is inconsistent with G.cells.indexMap.");

    result.cells.reserve(indexMap.size());
    std::unordered_set<CanonicalMeshIndex> cartesianIds;
    cartesianIds.reserve(indexMap.size());
    for (std::size_t storage = 0;
         storage < indexMap.size();
         ++storage)
    {
        if (indexMap[storage].empty())
            throw std::runtime_error(
                "G.cells.indexMap contains an empty row.");

        CanonicalMeshCell cell;
        cell.cartesianId =
            static_cast<CanonicalMeshIndex>(indexMap[storage][0]);
        if (!cartesianIds.insert(cell.cartesianId).second)
            throw std::runtime_error(
                "G.cells.indexMap contains duplicate Cartesian ids.");

        if (cellFacePos[storage].empty() ||
            cellFacePos[storage + 1].empty())
        {
            throw std::runtime_error(
                "G.cells.facePos contains an empty row.");
        }

        const int begin = cellFacePos[storage][0];
        const int end = cellFacePos[storage + 1][0];
        if (begin < 0 ||
            end < begin ||
            static_cast<std::size_t>(end) > cellFaces.size())
        {
            throw std::runtime_error(
                "Invalid cell-face offset range.");
        }

        cell.faces.reserve(
            static_cast<std::size_t>(end - begin));

        for (int entry = begin; entry < end; ++entry)
        {
            const auto &faceRow =
                cellFaces[static_cast<std::size_t>(entry)];
            if (faceRow.size() < 2)
                throw std::runtime_error(
                    "G.cells.faces row must contain face id and position.");

            const int inputFaceId = faceRow[0];
            const int inputPosition = faceRow[1];
            if (inputFaceId < 0 ||
                static_cast<std::size_t>(inputFaceId + 1) >= faceNodePos.size())
            {
                throw std::runtime_error(
                    "Cell references invalid face id.");
            }

            if (faceNodePos[static_cast<std::size_t>(inputFaceId)].empty() ||
                faceNodePos[static_cast<std::size_t>(inputFaceId + 1)].empty())
            {
                throw std::runtime_error(
                    "G.faces.nodePos contains an empty row.");
            }

            const int nodeBegin =
                faceNodePos[static_cast<std::size_t>(inputFaceId)][0];
            const int nodeEnd =
                faceNodePos[static_cast<std::size_t>(inputFaceId + 1)][0];
            if (nodeBegin < 0 ||
                nodeEnd < nodeBegin ||
                static_cast<std::size_t>(nodeEnd) > faceNodes.size())
            {
                throw std::runtime_error(
                    "Invalid face-node offset range.");
            }

            CanonicalMeshFace face;
            face.inputFaceId = static_cast<CanonicalMeshIndex>(inputFaceId);
            face.inputPosition = static_cast<CanonicalMeshIndex>(inputPosition);
            face.nodeIndices.reserve(
                static_cast<std::size_t>(nodeEnd - nodeBegin));

            for (int nodeEntry = nodeBegin;
                 nodeEntry < nodeEnd;
                 ++nodeEntry)
            {
                const auto &nodeRow =
                    faceNodes[static_cast<std::size_t>(nodeEntry)];
                if (nodeRow.empty())
                    throw std::runtime_error(
                        "G.faces.nodes contains an empty row.");

                const int nodeIndex = nodeRow[0];
                if (nodeIndex < 0 ||
                    static_cast<std::size_t>(nodeIndex) >= result.nodes.size())
                {
                    throw std::runtime_error(
                        "Face references invalid node id.");
                }

                face.nodeIndices.push_back(
                    static_cast<std::size_t>(nodeIndex));
            }

            cell.faces.push_back(std::move(face));
        }

        result.cells.push_back(std::move(cell));
    }

    if (result.nodes.empty() || result.cells.empty())
        throw std::runtime_error(
            "Mesh data contain no nodes or cells.");

    return result;
}

std::vector<std::array<CanonicalMeshIndex, 2>> loadMrstFaceNeighbors(
    const std::string &dataDirectory)
{
    const auto table =
        readMatlabIndexCsv<int>(
            path(dataDirectory, "G/faces/neighbors/data.csv"));

    std::vector<std::array<CanonicalMeshIndex, 2>> result;
    result.reserve(table.size());
    for (const auto &row : table)
    {
        if (row.size() < 2)
            throw std::runtime_error(
                "G.faces.neighbors contains an incomplete row.");
        result.push_back(
            {{static_cast<CanonicalMeshIndex>(row[0]),
              static_cast<CanonicalMeshIndex>(row[1])}});
    }
    return result;
}

} // namespace MPMC::detail
