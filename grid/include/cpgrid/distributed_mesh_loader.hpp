/**
 * @file distributed_mesh_loader.hpp
 * @brief CpGrid root-only 输入与 prepared Mesh 分发接口。
 */
#pragma once

#include <cpgrid/grdecl.hpp>
#include <cpgrid/mesh.hpp>

#include <memory>
#include <optional>
#include <string>

namespace MPMC
{

/** @brief root-only 网格输入的实际格式。 */
enum class DistributedMeshInputFormat
{
    MrstCsv,
    Grdecl
};

/**
 * @brief root-only ingest 的结果。
 *
 * `grdeclOnRoot` 只在 direct GRDECL 且当前 rank 为 root 时有值，用于后续
 * root-only 岩石属性装配和连接性报告；其它 rank 保持空值。
 */
struct DistributedMeshLoadResult final
{
    std::unique_ptr<Mesh> mesh;
    DistributedMeshInputFormat format{DistributedMeshInputFormat::MrstCsv};
    std::optional<GrdeclGridData> grdeclOnRoot;
};

/**
 * @brief 仅由 root 访问文件系统、解析输入并构造完整 Mesh，然后分发本地快照。
 *
 * `forceGrdecl=true` 对应显式 `-grdecl`。否则 root 根据输入路径是否为 `.grdecl`
 * 普通文件决定 GRDECL/MRST CSV 路径；判断结果会广播到全部 rank。
 */
[[nodiscard]] DistributedMeshLoadResult loadDistributedMeshFromRoot(
    const std::string &inputPath,
    bool forceGrdecl,
    const GrdeclLoadOptions &grdeclOptions = {},
    CpCommunicator comm = PETSC_COMM_WORLD,
    int root = 0);

} // namespace MPMC
