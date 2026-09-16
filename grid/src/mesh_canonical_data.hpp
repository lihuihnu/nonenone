/**
 * @file mesh_canonical_data.hpp
 * @brief Mesh 输入层使用的统一 canonical 中间表示。
 *
 * 该文件仅供 grid/src 内部使用。输入适配器负责把 MRST CSV 或 GRDECL
 * 转换为这里定义的稳定节点/单元/面索引表示；Mesh 本体随后只负责把该
 * 表示物化为 Node/Face/Polyhedron 拓扑对象。
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace MPMC
{

struct GrdeclGridData;

namespace detail
{

using CanonicalMeshIndex = std::int64_t;

struct CanonicalMeshNode final
{
    double x{0.0};
    double y{0.0};
    double z{0.0};
};

struct CanonicalMeshFace final
{
    std::vector<std::size_t> nodeIndices;
    CanonicalMeshIndex inputPosition{0};
    CanonicalMeshIndex inputFaceId{0};
};

struct CanonicalMeshCell final
{
    CanonicalMeshIndex cartesianId{0};
    std::vector<CanonicalMeshFace> faces;
};

struct CanonicalMeshData final
{
    std::string sourcePath;
    std::string sourceFormat;
    std::array<int, 3> logicalDimensions{{0, 0, 0}};

    std::vector<CanonicalMeshNode> nodes;
    std::vector<CanonicalMeshCell> cells;

    /** inputFaceId -> [first storage cell, second storage cell or -1]. */
    std::vector<std::array<CanonicalMeshIndex, 2>> faceNeighbors;
};

/** @brief 将 MRST CSV 主体数据转换为统一 canonical 网格表示。 */
[[nodiscard]] CanonicalMeshData loadMrstCanonicalMeshData(
    const std::string &dataDirectory);

/** @brief 延迟读取 MRST 邻接表，以保持 Mesh::initializeTopology() 语义。 */
[[nodiscard]] std::vector<std::array<CanonicalMeshIndex, 2>> loadMrstFaceNeighbors(
    const std::string &dataDirectory);

/** @brief 将已解析的 GRDECL 数据转换为统一 canonical 网格表示。 */
[[nodiscard]] CanonicalMeshData makeGrdeclCanonicalMeshData(
    const GrdeclGridData &grdecl);

} // namespace detail

} // namespace MPMC
