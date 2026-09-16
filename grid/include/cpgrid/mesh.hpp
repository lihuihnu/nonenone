/**
 * @file mesh.hpp
 * @brief CpGrid 网格拓扑、分区和原始输入编号的核心数据结构。
 */
#pragma once

#include <cpgrid/cartesian_cell_directory.hpp>
#include <cpgrid/config.hpp>
#include <cpgrid/filter_range.hpp>
#include <cpgrid/grdecl.hpp>
#include <cpgrid/node.hpp>
#include <cpgrid/polyhedron.hpp>
#include <cpgrid/types.hpp>

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace MPMC
{

namespace detail
{
struct CanonicalMeshData;
}

/**
 * @brief 非结构化/角点网格拓扑唯一拥有者。
 *
 * Mesh 拥有 Node 与 Polyhedron 的实际存储；Face 内的 Node* / Polyhedron*
 * 只观察这些对象。因此 Mesh 禁止复制和移动，且拓扑建立后禁止重新读取数据。
 */
class Mesh final
{
  public:
    explicit Mesh(
        std::string dataDirectory,
        CpCommunicator comm = PETSC_COMM_WORLD);

    /** @brief 由已解析的 GRDECL corner-point 数据直接构造 Mesh。 */
    explicit Mesh(
        const GrdeclGridData &grdecl,
        CpCommunicator comm = PETSC_COMM_WORLD);

    Mesh(const Mesh &) = delete;
    Mesh &operator=(const Mesh &) = delete;
    Mesh(Mesh &&) = delete;
    Mesh &operator=(Mesh &&) = delete;

    ~Mesh() = default;

    /**
     * @brief 集合式地把 root 上已完成 prepareForUse() 的完整 Mesh 分发为本地快照。
     *
     * root rank 必须传入非空且已完成 topology/partition/current-id 的 Mesh；
     * 其它 rank 必须传入空指针。函数保留 root 的完整 Mesh，并向其它 rank
     * 发送 owned+one-ring ghost 拓扑，同时广播轻量全局编号目录。
     */
    [[nodiscard]] static std::unique_ptr<Mesh> distributePreparedRoot(
        std::unique_ptr<Mesh> rootMesh,
        CpCommunicator comm = PETSC_COMM_WORLD,
        int root = 0);

    /**
     * @brief 建立 Face owner/neighbor 拓扑。
     *
     * 调用后 Node/Polyhedron 地址被视为永久稳定。
     */
    void initializeTopology();

    /**
     * @brief 使用 METIS 按 communicator rank 数进行分区。
     */
    void partition();

    /**
     * @brief 仅在指定 root 上执行 METIS，不调用任何分区广播 collective。
     *
     * 用于 root-only ingest：其它 rank 尚未构造完整 Mesh，因此不能参与历史
     * `partition()` 中的 owner broadcast。算法与 METIS 参数保持完全相同。
     */
    void partitionRootOnly(int root = 0);

    /**
     * @brief root-only ingest 的 topology -> partition -> current-id 准备路径。
     */
    void prepareForRootDistribution(int root = 0)
    {
        initializeTopology();
        partitionRootOnly(root);
        resetCurrentIds();
    }

    /**
     * @brief 按 rank 顺序重新分配连续 current cell id。
     *
     * 该编号用于 PETSc DOF 排列。
     */
    void resetCurrentIds();

    /**
     * @brief 完成 topology -> partition -> current-id 三阶段准备。
     */
    void prepareForUse()
    {
        initializeTopology();
        partition();
        resetCurrentIds();
    }

    /**
     * @brief 在非 root rank 上把 replicated 拓扑压缩为 owned+ghost 快照。
     *
     * 必须在 `prepareForUse()` 之后、创建 CpGrid/DofLayout 之前调用。rank 0
     * 默认保留完整拓扑，用于全局编号表和可选 VTK 输出；其它 rank 释放远端
     * Node/Face/Polyhedron，仅保留求解所需的一环局部拓扑。全局编号目录保留。
     */
    void compactToLocalSnapshot(bool keepFullMeshOnRoot = true);

    [[nodiscard]] bool isLocalSnapshotCompacted() const noexcept
    {
        return localSnapshotCompacted_;
    }

    [[nodiscard]] bool
    isTopologyInitialized() const noexcept
    {
        return topologyInitialized_;
    }

    [[nodiscard]] bool
    isReadyForDofs() const noexcept
    {
        return topologyInitialized_ &&
               partitioned_ &&
               currentIdsReady_;
    }

    [[nodiscard]] CpCommunicator
    communicator() const noexcept
    {
        return comm_;
    }

    [[nodiscard]] int processCount() const noexcept;
    [[nodiscard]] int rank() const noexcept;

    [[nodiscard]] std::size_t
    nodeCount() const noexcept
    {
        return globalNodeCount_;
    }

    [[nodiscard]] std::size_t
    cellCount() const noexcept
    {
        return globalCellCount_;
    }

    /** @brief 当前 rank 实际物化的拓扑单元数；replicated 模式等于 cellCount()。 */
    [[nodiscard]] std::size_t materializedCellCount() const noexcept
    {
        return cells_.size();
    }

    /** @brief 当前 rank 实际物化的节点数；replicated 模式等于 nodeCount()。 */
    [[nodiscard]] std::size_t materializedNodeCount() const noexcept
    {
        return nodes_.size();
    }

    [[nodiscard]] PetscInt
    localCellCount() const noexcept
    {
        return localCellCount_;
    }

    /** @brief 当前 rank 在 rank-major current-id 中的首个 owned cell id。 */
    [[nodiscard]] PetscInt firstOwnedCurrentId() const
    {
        requireCurrentIds_();
        return rankCurrentOffsets_.at(static_cast<std::size_t>(rank()));
    }

    /** @brief 仅凭全局 current id 查询 owner rank，不要求远端 topology 已物化。 */
    [[nodiscard]] int ownerRankFromCurrentId(PetscInt currentId) const
    {
        requireCurrentIds_();
        if (currentId < 0 || static_cast<std::size_t>(currentId) >= globalCellCount_)
            throw std::out_of_range("Current cell id is out of range.");
        const auto upper = std::upper_bound(
            rankCurrentOffsets_.begin(),
            rankCurrentOffsets_.end(),
            currentId);
        if (upper == rankCurrentOffsets_.begin() || upper == rankCurrentOffsets_.end())
            throw std::logic_error("Mesh rank current-id offsets are inconsistent.");
        return static_cast<int>((upper - rankCurrentOffsets_.begin()) - 1);
    }

    /**
     * @brief 当前 rank 拥有的 current cell ids。
     *
     * resetCurrentIds() 后按 current-id 升序缓存；该顺序等价于历史
     * `localCells()` 的 rank 内稳定 storage 顺序。
     */
    [[nodiscard]] const std::vector<PetscInt> &
    ownedCellIds() const
    {
        requireCurrentIds_();
        return ownedCellIds_;
    }

    /** @brief 当前 rank 求解所需的一环 ghost current cell ids。 */
    [[nodiscard]] const std::vector<PetscInt> &
    ghostCellIds() const
    {
        requireCurrentIds_();
        return ghostCellIds_;
    }

    [[nodiscard]] const std::string &sourcePath() const noexcept
    {
        return sourcePath_;
    }

    [[nodiscard]] const std::string &sourceFormat() const noexcept
    {
        return sourceFormat_;
    }

    [[nodiscard]] std::array<int, 3> logicalDimensions() const noexcept
    {
        return logicalDimensions_;
    }

    [[nodiscard]] std::size_t faceInstanceCount() const noexcept
    {
        return globalFaceInstanceCount_;
    }

    [[nodiscard]] std::size_t uniqueFaceCount() const noexcept
    {
        return globalUniqueFaceCount_ > 0 ?
            globalUniqueFaceCount_ :
            globalFaceInstanceCount_;
    }

    /**
     * @brief 按原始网格数据文件中的行序号访问单元。
     *
     * The input index is the zero-based row number used by
     * `G/cells/indexMap/data.csv` and by the rock CSV files.  It never changes
     * after METIS partitioning and is the canonical id exposed to cases, wells
     * and post-processing.
     */
    [[nodiscard]] Polyhedron &cellByInputIndex(std::size_t index);

    [[nodiscard]] const Polyhedron &cellByInputIndex(std::size_t index) const;

    /** @brief 将原始输入单元序号转换为分区后的 current id。 */
    [[nodiscard]] PetscInt currentIdFromInputIndex(std::size_t index) const
    {
        requireCurrentIds_();
        return currentIdByInputIndex_.at(index);
    }

    /** @brief 将分区后的 current id 转回原始输入单元序号。 */
    [[nodiscard]] PetscInt inputIndexFromCurrentId(PetscInt currentId) const
    {
        requireCurrentIds_();
        if (currentId < 0 ||
            static_cast<std::size_t>(currentId) >= inputIndexByCurrentId_.size())
        {
            throw std::out_of_range("Current cell id is out of range.");
        }
        return inputIndexByCurrentId_[static_cast<std::size_t>(currentId)];
    }

    /** @brief 将 0 基 Cartesian id 直接解析为 current id，不要求远端 cell 已物化。 */
    [[nodiscard]] PetscInt currentIdFromCartesianId(PetscInt cartesianId) const
    {
        requireCurrentIds_();
        return currentIdFromInputIndex(
            static_cast<std::size_t>(cartesianDirectory_.inputIndex(cartesianId)));
    }

    /**
     * @brief 按内部存储序号访问单元；数值上与 input index 一致。
     *
     * This alias remains for mesh implementation code.  User-facing code should
     * use `cellByInputIndex()` to make the numbering convention explicit.
     */
    [[nodiscard]] Polyhedron &
    cellByStorageIndex(std::size_t index)
    {
        return cellByInputIndex(index);
    }

    [[nodiscard]] const Polyhedron &
    cellByStorageIndex(std::size_t index) const
    {
        return cellByInputIndex(index);
    }

    /**
     * @brief 按 resetCurrentIds() 后的 current id 访问单元。
     */
    [[nodiscard]] Polyhedron &
    cellByCurrentId(PetscInt currentId);

    [[nodiscard]] const Polyhedron &
    cellByCurrentId(PetscInt currentId) const;

    /**
     * @brief 按 MRST G.cells.indexMap 对应的 0 基 Cartesian id 访问单元。
     */
    [[nodiscard]] Polyhedron &
    cellByCartesianId(PetscInt cartesianId);

    [[nodiscard]] const Polyhedron &
    cellByCartesianId(PetscInt cartesianId) const;

    [[nodiscard]] const Node &node(std::size_t index) const;

    [[nodiscard]] auto allCells() noexcept
    {
        return FilterRange(
            cells_.begin(),
            cells_.end(),
            [](const Polyhedron &)
            {
                return true;
            });
    }

    [[nodiscard]] auto allCells() const noexcept
    {
        return FilterRange(
            cells_.cbegin(),
            cells_.cend(),
            [](const Polyhedron &)
            {
                return true;
            });
    }

    [[nodiscard]] auto localCells()
    {
        requirePartitioned_();
        const int ownerRank = rank();

        return FilterRange(
            cells_.begin(),
            cells_.end(),
            [ownerRank](const Polyhedron &cell)
            {
                return cell.isOwnedBy(ownerRank);
            });
    }

    [[nodiscard]] auto localCells() const
    {
        requirePartitioned_();
        const int ownerRank = rank();

        return FilterRange(
            cells_.cbegin(),
            cells_.cend(),
            [ownerRank](const Polyhedron &cell)
            {
                return cell.isOwnedBy(ownerRank);
            });
    }

    /**
     * @brief 输出基础网格摘要。
     */

    /**
     * @brief 保存 input/cartesian/current 单元编号映射用于诊断。
     */
    void writeCellIdMap(
        const std::string &filename) const;

  private:
    /** @brief 仅供 distributePreparedRoot() 在非 root 构造空接收壳。 */
    explicit Mesh(CpCommunicator comm);

    void initializeCommunicator_();
    void initializeFromCanonical_(const detail::CanonicalMeshData &data);
    void bindNeighbors_();
    void rebuildCurrentIdLookup_();
    void rebuildLocalPartitionView_();
    void validateMaterializedOrdering_() const;
    [[nodiscard]] std::size_t materializedStorageIndexFromCurrentId_(
        PetscInt currentId) const;

    void requireTopology_() const;
    void requirePartitioned_() const;
    void requireCurrentIds_() const;

    [[nodiscard]] std::string
    path_(const std::string &relative) const;

    CpCommunicator comm_{PETSC_COMM_WORLD};
    int processCount_{0};
    int rank_{0};
    std::string dataDirectory_;
    std::string sourcePath_;
    std::string sourceFormat_{"MRST_CSV"};
    std::array<int, 3> logicalDimensions_{{0, 0, 0}};

    std::size_t globalNodeCount_{0};
    std::size_t globalCellCount_{0};
    std::size_t globalFaceInstanceCount_{0};
    std::size_t globalUniqueFaceCount_{0};

    /** inputFaceId -> [first storage cell, second storage cell or -1]. */
    std::vector<std::array<PetscInt, 2>> faceNeighbors_;

    std::vector<Node> nodes_;
    std::vector<Polyhedron> cells_;

    /** 全局 Cartesian id -> input index 的紧凑只读目录。 */
    detail::CompactCartesianCellDirectory<PetscInt> cartesianDirectory_;

    /** 全局 input index -> current id 的轻量编号目录。 */
    std::vector<PetscInt> currentIdByInputIndex_;

    /** 全局 current id -> input index 的轻量编号目录。 */
    std::vector<PetscInt> inputIndexByCurrentId_;

    /** rank-major current-id 区间前缀；长度 processCount+1。 */
    std::vector<PetscInt> rankCurrentOffsets_;

    /** 当前 rank 的 owned current ids，按 current id 升序。 */
    std::vector<PetscInt> ownedCellIds_;

    /** 当前 rank 的一环 ghost current ids，升序且唯一。 */
    std::vector<PetscInt> ghostCellIds_;

    /** compact 模式下验证 owned/ghost/node 排序不变量；不再保存持久 hash 表。 */

    PetscInt localCellCount_{0};
    bool localSnapshotCompacted_{false};

    bool topologyInitialized_{false};
    bool partitioned_{false};
    bool currentIdsReady_{false};
};

} // namespace MPMC
