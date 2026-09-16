/**
 * @file cpgrid.hpp
 * @brief grid 模块中的 `cpgrid` 源码。
 */
#pragma once

#include <cpgrid/layout_registry.hpp>
#include <cpgrid/mesh.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief 非结构化/角点网格的 PETSc 离散接口。
 *
 * CpGrid 只负责：
 * - Mesh 拓扑访问；
 * - per-cell DOF 布局；
 * - PETSc Vec/Mat/DM 创建；
 * - 孔隙度/渗透率场；
 * - 几何量与 TPFA 几何传递系数缓存。
 *
 * CpGrid 不知道相、组分、EOS、flash、SCW 或化学反应。
 *
 * 每单元 DOF 布局由上层模型显式调用 `registerLayout()` 注册；CpGridCore
 * 自身只预注册 permeability(3) 与 porosity(1)，不理解任何物理/AD 语义。
 */
class CpGridCore
{
  public:
    using Cell = Polyhedron;
    using FaceType = Face;

    /**
     * @brief 只读本地岩石属性快照。
     *
     * 构造时同步 porosity/permeability 的 owned + ghost 数据；
     * 局部 Vec 从 PETSc DM 临时向量池借用，析构时自动恢复数组并归还。
     */
    class RockLocalView final
    {
      public:
        explicit RockLocalView(
            const CpGridCore &grid)
            : grid_(grid)
        {
            permeabilityLocal_ =
                grid_.borrowLocalVector(
                    3,
                    grid_.permeability_);

            porosityLocal_ =
                grid_.borrowLocalVector(
                    1,
                    grid_.porosity_);

            PetscCallAbort(
                PETSC_COMM_SELF,
                VecGetArrayRead(
                    permeabilityLocal_,
                    &permeabilityArray_));

            PetscCallAbort(
                PETSC_COMM_SELF,
                VecGetArrayRead(
                    porosityLocal_,
                    &porosityArray_));
        }

        RockLocalView(
            const RockLocalView &) = delete;
        RockLocalView &operator=(
            const RockLocalView &) = delete;

        RockLocalView(
            RockLocalView &&other) noexcept
            : grid_(other.grid_),
              permeabilityLocal_(
                  std::exchange(
                      other.permeabilityLocal_,
                      nullptr)),
              porosityLocal_(
                  std::exchange(
                      other.porosityLocal_,
                      nullptr)),
              permeabilityArray_(
                  std::exchange(
                      other.permeabilityArray_,
                      nullptr)),
              porosityArray_(
                  std::exchange(
                      other.porosityArray_,
                      nullptr))
        {
        }

        RockLocalView &operator=(
            RockLocalView &&) = delete;

        ~RockLocalView() noexcept
        {
            release_();
        }

        /**
         * @brief 读取本地或 ghost 单元孔隙度。
         */
        [[nodiscard]] double porosity(
            const Cell &cell) const
        {
            requireValid_();

            const PetscInt localIndex =
                grid_.dofMap(1)
                    .localIndex(
                        cell,
                        0);

            const double value =
                PetscRealPart(
                    porosityArray_[
                        localIndex]);

            if (!std::isfinite(value) ||
                value < 0.0 ||
                value > 1.0)
            {
                throw std::runtime_error(
                    "CpGrid porosity must be finite and inside [0,1].");
            }

            return value;
        }

        /**
         * @brief 读取本地或 ghost 单元的三个主方向渗透率。
         */
        [[nodiscard]] std::array<double, 3>
        permeability(
            const Cell &cell) const
        {
            requireValid_();

            const DofMap &map =
                grid_.dofMap(3);

            std::array<double, 3> value{};

            for (PetscInt axis = 0;
                 axis < 3;
                 ++axis)
            {
                const PetscInt localIndex =
                    map.localIndex(
                        cell,
                        axis);

                value[
                    static_cast<std::size_t>(
                        axis)] =
                    PetscRealPart(
                        permeabilityArray_[
                            localIndex]);
            }

            return value;
        }

      private:
        void requireValid_() const
        {
            if (permeabilityLocal_ == nullptr ||
                porosityLocal_ == nullptr ||
                permeabilityArray_ == nullptr ||
                porosityArray_ == nullptr)
            {
                throw std::logic_error(
                    "RockLocalView has already released its PETSc resources.");
            }
        }

        void release_() noexcept
        {
            if (permeabilityLocal_ != nullptr &&
                permeabilityArray_ != nullptr)
            {
                PetscCallAbort(
                    PETSC_COMM_SELF,
                    VecRestoreArrayRead(
                        permeabilityLocal_,
                        &permeabilityArray_));
            }

            if (porosityLocal_ != nullptr &&
                porosityArray_ != nullptr)
            {
                PetscCallAbort(
                    PETSC_COMM_SELF,
                    VecRestoreArrayRead(
                        porosityLocal_,
                        &porosityArray_));
            }

            if (permeabilityLocal_ != nullptr)
                grid_.restoreLocalVector(3, permeabilityLocal_);

            if (porosityLocal_ != nullptr)
                grid_.restoreLocalVector(1, porosityLocal_);

            permeabilityArray_ = nullptr;
            porosityArray_ = nullptr;
        }

        const CpGridCore &grid_;

        Vec permeabilityLocal_{nullptr};
        Vec porosityLocal_{nullptr};

        const PetscScalar *
            permeabilityArray_{nullptr};
        const PetscScalar *
            porosityArray_{nullptr};
    };

    /**
     * @brief 从已经 prepareForUse() 的 Mesh 创建 CpGrid。
     *
     * Mesh 生命周期必须覆盖 CpGrid。
     */
    explicit CpGridCore(
        Mesh &mesh,
        PetscInt legacyPrimaryDof = 0)
        : mesh_(mesh),
          layouts_(mesh)
    {
        if (!mesh_.isReadyForDofs())
        {
            throw std::invalid_argument(
                "CpGrid requires Mesh::prepareForUse() first.");
        }

        registerStandardLayouts_(legacyPrimaryDof);

        permeability_ =
            createGlobalVector(3);
        porosity_ =
            createGlobalVector(1);
    }

    CpGridCore(const CpGridCore &) = delete;
    CpGridCore &operator=(
        const CpGridCore &) = delete;
    CpGridCore(CpGridCore &&) = delete;
    CpGridCore &operator=(
        CpGridCore &&) = delete;

    ~CpGridCore() noexcept
    {
        releaseRockVectors_();
    }

    /** @brief 访问外部拥有的 Mesh 拓扑与分区；CpGrid 不接管其生命周期。 */
    [[nodiscard]] Mesh &mesh() noexcept
    {
        return mesh_;
    }

    /** @brief 只读访问底层 Mesh 拓扑和 communicator。 */
    [[nodiscard]] const Mesh &mesh() const noexcept
    {
        return mesh_;
    }

    /**
     * @brief 显式注册额外的每单元 DOF 布局。
     *
     * 该操作包含 MPI collective，所有 rank 必须以一致顺序调用。
     */
    DofLayout &registerLayout(
        PetscInt dofPerCell)
    {
        return layouts_.registerLayout(
            dofPerCell);
    }

    /** @brief 查询指定每单元 DOF 布局是否已注册。 */
    [[nodiscard]] bool hasLayout(
        PetscInt dofPerCell) const noexcept
    {
        return layouts_.contains(
            dofPerCell);
    }

    /** @brief 按每单元 DOF 数访问已注册的 PETSc 布局。 */
    [[nodiscard]] DofLayout &layout(
        PetscInt dofPerCell)
    {
        return layouts_.layout(
            dofPerCell);
    }

    /** @brief 只读访问已注册的每单元 PETSc 布局。 */
    [[nodiscard]] const DofLayout &layout(
        PetscInt dofPerCell) const
    {
        return layouts_.layout(
            dofPerCell);
    }

    /** @brief 访问已注册布局对应的 owned/ghost DOF 映射。 */
    [[nodiscard]] DofMap &dofMap(
        PetscInt dofPerCell)
    {
        return layout(
                   dofPerCell)
            .dofMap();
    }

    /** @brief 只读访问 owned/ghost DOF 映射。 */
    [[nodiscard]] const DofMap &dofMap(
        PetscInt dofPerCell) const
    {
        return layout(
                   dofPerCell)
            .dofMap();
    }

    /** @brief 借用表示已注册布局的 PETSc DM；调用方不得销毁。 */
    [[nodiscard]] DM dm(
        PetscInt dofPerCell) const
    {
        return layout(
                   dofPerCell)
            .dm();
    }

    /** @brief 创建每个 owned 单元一个连续 block 的 PETSc 全局向量。 */
    [[nodiscard]] Vec createGlobalVector(
        PetscInt dofPerCell) const
    {
        return layout(
                   dofPerCell)
            .createGlobalVector();
    }

    /**
     * @brief 创建按原始输入文件单元顺序排列的全局向量副本。
     *
     * CpGrid solves internally with `current id` ordering because PETSc requires each
     * MPI rank to own a contiguous global range.  `Mesh::resetCurrentIds()` therefore
     * reorders cells by MPI owner.  That internal ordering is intentionally kept for
     * assembly and communication.
     *
     * This routine performs the inverse permutation only for external data:
     *
     * ```text
     * internal block = current_id(cell)
     * output   block = input_index(cell)
     * ```
     *
     * Each block still stores its `dofPerCell` components contiguously, so the exported
     * vector layout is exactly:
     *
     * ```text
     * input cell 0: component 0 ... component dof-1
     * input cell 1: component 0 ... component dof-1
     * ...
     * ```
     *
     * The returned Vec is newly allocated and must be destroyed by the caller.
     */
    [[nodiscard]] Vec createInputOrderedCopy(
        PetscInt dofPerCell,
        Vec currentOrderedVector) const
    {
        if (currentOrderedVector == nullptr)
        {
            throw std::invalid_argument(
                "createInputOrderedCopy requires a valid PETSc Vec.");
        }

        const DofMap &map = dofMap(dofPerCell);

        PetscInt globalSize = 0;
        PetscCallAbort(
            mesh_.communicator(),
            VecGetSize(currentOrderedVector, &globalSize));

        if (globalSize != map.globalDofCount())
        {
            throw std::invalid_argument(
                "Vector size does not match the requested CpGrid DOF layout.");
        }

        PetscInt localSize = 0;
        PetscCallAbort(
            mesh_.communicator(),
            VecGetLocalSize(currentOrderedVector, &localSize));

        Vec inputOrdered = nullptr;
        PetscCallAbort(
            mesh_.communicator(),
            VecCreateMPI(
                mesh_.communicator(),
                localSize,
                globalSize,
                &inputOrdered));
        PetscCallAbort(
            mesh_.communicator(),
            VecSetBlockSize(inputOrdered, dofPerCell));
        PetscCallAbort(
            mesh_.communicator(),
            VecSet(inputOrdered, 0.0));

        PetscInt sourceBegin = 0;
        PetscInt sourceEnd = 0;
        const PetscScalar *source = nullptr;
        PetscCallAbort(
            mesh_.communicator(),
            VecGetOwnershipRange(
                currentOrderedVector,
                &sourceBegin,
                &sourceEnd));
        PetscCallAbort(
            mesh_.communicator(),
            VecGetArrayRead(
                currentOrderedVector,
                &source));

        std::vector<PetscInt> destinationIndices(
            static_cast<std::size_t>(dofPerCell));
        std::vector<PetscScalar> values(
            static_cast<std::size_t>(dofPerCell));

        for (PetscInt currentId : mesh_.ownedCellIds())
        {
            const Cell &cell = mesh_.cellByCurrentId(currentId);
            const PetscInt inputCell = cell.inputIndex();

            for (PetscInt component = 0;
                 component < dofPerCell;
                 ++component)
            {
                const PetscInt sourceGlobal =
                    map.globalIndex(cell, component);

                if (sourceGlobal < sourceBegin ||
                    sourceGlobal >= sourceEnd)
                {
                    PetscCallAbort(
                        mesh_.communicator(),
                        VecRestoreArrayRead(
                            currentOrderedVector,
                            &source));
                    PetscCallAbort(
                        mesh_.communicator(),
                        VecDestroy(&inputOrdered));
                    throw std::runtime_error(
                        "Owned CpGrid cell is outside the source Vec ownership range.");
                }

                destinationIndices[
                    static_cast<std::size_t>(component)] =
                    inputCell * dofPerCell + component;

                values[
                    static_cast<std::size_t>(component)] =
                    source[sourceGlobal - sourceBegin];
            }

            PetscCallAbort(
                mesh_.communicator(),
                VecSetValues(
                    inputOrdered,
                    dofPerCell,
                    destinationIndices.data(),
                    values.data(),
                    INSERT_VALUES));
        }

        PetscCallAbort(
            mesh_.communicator(),
            VecRestoreArrayRead(
                currentOrderedVector,
                &source));
        PetscCallAbort(
            mesh_.communicator(),
            VecAssemblyBegin(inputOrdered));
        PetscCallAbort(
            mesh_.communicator(),
            VecAssemblyEnd(inputOrdered));

        return inputOrdered;
    }

    /** @brief 为已注册布局创建 owned+ghost 局部向量。 */
    [[nodiscard]] Vec createLocalVector(
        PetscInt dofPerCell) const
    {
        return layout(
                   dofPerCell)
            .createLocalVector();
    }

    /** @brief 创建局部向量并立即把 `global` scatter 到 owned+ghost 布局。 */
    [[nodiscard]] Vec createLocalVector(
        PetscInt dofPerCell,
        Vec global) const
    {
        return layout(
                   dofPerCell)
            .createLocalVector(
                global);
    }

    /** @brief 从 PETSc DM 临时向量池借用 owned+ghost 局部向量。 */
    [[nodiscard]] Vec borrowLocalVector(
        PetscInt dofPerCell) const
    {
        return layout(dofPerCell).borrowLocalVector();
    }

    /** @brief 借用局部向量并同步 global 数据。 */
    [[nodiscard]] Vec borrowLocalVector(
        PetscInt dofPerCell,
        Vec global) const
    {
        return layout(dofPerCell).borrowLocalVector(global);
    }

    /** @brief 将临时局部向量归还对应布局的 PETSc DM 池。 */
    void restoreLocalVector(
        PetscInt dofPerCell,
        Vec &local) const noexcept
    {
        layout(dofPerCell).restoreLocalVector(local);
    }

    /** @brief 将全局单元向量 scatter 到 owned+ghost 局部布局。 */
    void globalToLocal(
        PetscInt dofPerCell,
        Vec global,
        Vec local,
        InsertMode mode = INSERT_VALUES) const
    {
        layout(
            dofPerCell)
            .globalToLocal(
                global,
                local,
                mode);
    }

    /** @brief 将局部 owned+ghost 向量累积回全局布局。 */
    void localToGlobal(
        PetscInt dofPerCell,
        Vec local,
        Vec global,
        InsertMode mode = ADD_VALUES) const
    {
        layout(
            dofPerCell)
            .localToGlobal(
                local,
                global,
                mode);
    }

    /** @brief 根据网格邻接关系预分配并创建 block sparse PETSc 矩阵。 */
    [[nodiscard]] Mat createMatrix(
        PetscInt dofPerCell) const
    {
        return layout(
                   dofPerCell)
            .createMatrix();
    }

    /**
     * @brief 当前 rank 拥有的单元范围。
     */
    [[nodiscard]] auto localCells()
    {
        return mesh_.localCells();
    }

    [[nodiscard]] auto localCells() const
    {
        return mesh_.localCells();
    }

    /** @brief 访问与指定网格单元相连的面。 */
    [[nodiscard]] std::vector<Face> &
    faces(Cell &cell) noexcept
    {
        return cell.faces();
    }

    /** @brief 只读访问与指定网格单元相连的面。 */
    [[nodiscard]] const std::vector<Face> &
    faces(const Cell &cell) const noexcept
    {
        return cell.faces();
    }

    /**
     * @brief 返回内部面的邻接单元。
     *
     * @throws std::invalid_argument 边界面没有邻居时抛出。
     */
    [[nodiscard]] Cell &neighbor(
        FaceType &face)
    {
        Cell *cell =
            face.neighborCell();

        if (cell == nullptr)
        {
            throw std::invalid_argument(
                "Boundary face has no neighbor cell.");
        }

        return *cell;
    }

    [[nodiscard]] const Cell &neighbor(
        const FaceType &face) const
    {
        const Cell *cell =
            face.neighborCell();

        if (cell == nullptr)
        {
            throw std::invalid_argument(
                "Boundary face has no neighbor cell.");
        }

        return *cell;
    }

    /**
     * @brief 返回已缓存的内部面 TPFA 几何传递系数。
     *
     * The cache contains the harmonic combination of the two half-face terms,
     * `T = (1/T_i + 1/T_j)^{-1}`; phase mobility is applied by the flow model.
     */
    [[nodiscard]] double transmissibility(
        const FaceType &face) const noexcept
    {
        return face.transmissibility();
    }

    /** @brief 返回单元几何 bulk volume [m^3]。 */
    [[nodiscard]] double cellVolume(
        const Cell &cell) const noexcept
    {
        return cell.volume();
    }

    /** @brief 返回面两侧已缓存的重力势差 `g (z_i-z_j)`。 */
    [[nodiscard]] double gravityTerm(
        const FaceType &face) const noexcept
    {
        return face.gravityTerm();
    }

    /** @brief 判断该单元是否由当前 MPI rank 拥有而非仅为 ghost。 */
    [[nodiscard]] bool isLocal(
        const Cell &cell) const noexcept
    {
        return cell.isOwnedBy(
            mesh_.rank());
    }

    /**
     * @brief 获取 CpGrid 拥有的全局渗透率 Vec。
     *
     * 布局为每单元 3 DOF：Kx、Ky、Kz。
     * CpGrid 不定义外部数据单位；求解器项目应统一采用 SI。
     */
    [[nodiscard]] Vec permeabilityVector() const noexcept
    {
        return permeability_;
    }

    /**
     * @brief 获取 CpGrid 拥有的全局孔隙度 Vec。
     */
    [[nodiscard]] Vec porosityVector() const noexcept
    {
        return porosity_;
    }

    /**
     * @brief 创建当前岩石属性的只读 owned+ghost 快照。
     */
    [[nodiscard]] RockLocalView
    rockLocalView() const
    {
        return RockLocalView(*this);
    }

    /**
     * @brief 根据当前 Mesh 与 permeability Vec 建立几何缓存。
     *
     * 每个 rank 只计算 owned + one-ring ghost 的 cell/face 纯几何；
     * 只有 owned cell faces 计算 transmissibility / gravity term。
     */
    void setup()
    {
        if (setupDone_)
        {
            throw std::logic_error(
                "CpGrid::setup may only be called once.");
        }

        cachePureGeometry_();

        const RockLocalView rock =
            rockLocalView();

        cacheLocalFlowGeometry_(
            rock);

        setupDone_ = true;
    }

    /** @brief 判断 `setup()` 是否已构建几何与流动缓存。 */
    [[nodiscard]] bool isSetup() const noexcept
    {
        return setupDone_;
    }

  protected:
    /** @brief 为单个已物化单元缓存纯几何；兼容 façade 可复用该入口。 */
    void cacheCellPureGeometry_(Cell &cell)
    {
        const PolyhedronGeometry cellGeometry =
            cell.computeGeometry();
        const double volume = cellGeometry.volume;
        const Point &centroid = cellGeometry.centroid;

        if (!std::isfinite(volume) ||
            volume <= 0.0 ||
            !centroid.isFinite())
        {
            throw std::runtime_error(
                "CpGrid encountered invalid cell geometry.");
        }

        cell.cacheVolume(volume);
        cell.cacheCentroid(centroid);

        for (FaceType &face : cell.faces())
        {
            const PolygonGeometry faceGeometry =
                face.computeGeometry();
            const double area = faceGeometry.area;
            const Point &faceCentroid = faceGeometry.centroid;
            const Point &unitNormal = faceGeometry.unitNormal;

            if (!std::isfinite(area) ||
                area <= 0.0 ||
                !faceCentroid.isFinite() ||
                !unitNormal.isFinite())
            {
                throw std::runtime_error(
                    "CpGrid encountered invalid face geometry.");
            }

            face.cacheArea(area);
            face.cacheCentroid(faceCentroid);
            face.cacheUnitNormal(unitNormal);

            /*
             * 当前基础模型对外边界采用 no-flow 几何传递系数。
             * 需要其它边界条件时由残差边界项显式处理。
             */
            if (face.isBoundary())
            {
                face.cacheTransmissibility(0.0);
                face.cacheGravityTerm(0.0);
            }
        }
    }

  private:
    void registerStandardLayouts_(PetscInt legacyPrimaryDof)
    {
        // 新核心仅需要 permeability(3) 与 porosity(1)。兼容 façade 可通过
        // legacyPrimaryDof 在它们之前注册旧主变量布局，以保持历史 collective 顺序。
        if (legacyPrimaryDof > 0)
            registerLayout(legacyPrimaryDof);

        registerLayout(3);
        registerLayout(1);
    }

    void cachePureGeometry_()
    {
        /*
         * 流动装配只需要 owned cells 及其一环 ghost。当前仍保留 replicated
         * Mesh 作为兼容/输出后备，但生产 Core 不再在每个 rank 上预计算
         * 完整网格的纯几何。
         */
        for (PetscInt cellId : mesh_.ownedCellIds())
            cacheCellPureGeometry_(mesh_.cellByCurrentId(cellId));
        for (PetscInt cellId : mesh_.ghostCellIds())
            cacheCellPureGeometry_(mesh_.cellByCurrentId(cellId));
    }

    void cacheLocalFlowGeometry_(
        const RockLocalView &rock)
    {
        for (PetscInt cellId : mesh_.ownedCellIds())
        {
            Cell &cell = mesh_.cellByCurrentId(cellId);
            for (FaceType &face :
                 cell.faces())
            {
                if (face.isBoundary())
                {
                    continue;
                }

                const Cell &adjacent =
                    neighbor(face);

                face.cacheTransmissibility(
                    computeTransmissibility_(
                        cell,
                        adjacent,
                        face,
                        rock));

                face.cacheGravityTerm(
                    computeGravityTerm_(
                        cell,
                        adjacent));
            }
        }
    }

    /**
     * @brief 计算从单元中心到面中心的 TPFA 半面传递系数。
     *
     * For diagonal permeability `K`, face area `A`, unit normal `n` and
     * centre-to-face vector `d`, this implementation uses
     * `T_half = A |d . (K n)| / |d|^2`. Two neighboring half terms are later
     * combined harmonically. Permeability is SI [m^2], so `T` has units [m^3].
     */
    [[nodiscard]] double
    computeHalfTransmissibility_(
        const Cell &cell,
        const FaceType &face,
        const RockLocalView &rock) const
    {
        const Point &cellCenter =
            cell.cachedCentroid();

        const Point &faceCenter =
            face.cachedCentroid();

        const Point distance =
            faceCenter -
            cellCenter;

        const double distanceSquared =
            distance.normSquared();

        if (!std::isfinite(
                distanceSquared) ||
            distanceSquared <= 0.0)
        {
            throw std::runtime_error(
                "Cell center and face center have zero or invalid distance.");
        }

        const Point &normal =
            face.cachedUnitNormal();

        const auto permeability =
            rock.permeability(cell);

        for (double value :
             permeability)
        {
            if (!std::isfinite(value) ||
                value < 0.0)
            {
                throw std::runtime_error(
                    "CpGrid permeability must be finite and non-negative.");
            }
        }

        const Point kNormal(
            permeability[0] *
                normal.x(),
            permeability[1] *
                normal.y(),
            permeability[2] *
                normal.z());

        const double halfTrans =
            face.area() *
            std::abs(
                distance.dot(kNormal)) /
            distanceSquared;

        if (!std::isfinite(halfTrans) ||
            halfTrans < 0.0)
        {
            throw std::runtime_error(
                "CpGrid produced an invalid half transmissibility.");
        }

        return halfTrans;
    }

    [[nodiscard]] double
    computeTransmissibility_(
        const Cell &cell,
        const Cell &adjacent,
        const FaceType &face,
        const RockLocalView &rock) const
    {
        const double first =
            computeHalfTransmissibility_(
                cell,
                face,
                rock);

        const FaceType *adjacentFace =
            face.neighborFace();
        if (adjacentFace == nullptr)
        {
            throw std::runtime_error(
                "Internal CpGrid face is missing its opposite face geometry.");
        }

        const double second =
            computeHalfTransmissibility_(
                adjacent,
                *adjacentFace,
                rock);

        if (first <= 0.0 ||
            second <= 0.0)
        {
            return 0.0;
        }

        const double denominator =
            1.0 / first +
            1.0 / second;

        if (!std::isfinite(denominator) ||
            denominator <= 0.0)
        {
            throw std::runtime_error(
                "CpGrid transmissibility denominator is invalid.");
        }

        const double transmissibility =
            1.0 / denominator;

        if (!std::isfinite(
                transmissibility) ||
            transmissibility < 0.0)
        {
            throw std::runtime_error(
                "CpGrid produced an invalid face transmissibility.");
        }

        return transmissibility;
    }

    [[nodiscard]] static double
    computeGravityTerm_(
        const Cell &cell,
        const Cell &adjacent)
    {
        constexpr double gravity =
            9.80665;

        return (
                   cell.cachedCentroid().z() -
                   adjacent.cachedCentroid().z()) *
               gravity;
    }

    void releaseRockVectors_() noexcept
    {
        if (permeability_ != nullptr)
        {
            PetscCallAbort(
                mesh_.communicator(),
                VecDestroy(
                    &permeability_));
        }

        if (porosity_ != nullptr)
        {
            PetscCallAbort(
                mesh_.communicator(),
                VecDestroy(
                    &porosity_));
        }
    }

    Mesh &mesh_;
    LayoutRegistry layouts_;

    Vec permeability_{nullptr};
    Vec porosity_{nullptr};

    bool setupDone_{false};
 };

/**
 * @brief 兼容旧 `CpGrid<Tag>` 源码接口的薄 façade。
 *
 * `CpGridCore` 不依赖任何模型 Tag；本层只在构造阶段把旧主变量 DOF 传给
 * Core，以保持历史布局注册顺序。新代码应使用 `CpGridCore` 并由模型层显式
 * 注册主变量/phase-state 等布局。
 */
template <class Tag>
class CpGrid final : public CpGridCore
{
  public:
    explicit CpGrid(Mesh &mesh)
        : CpGridCore(mesh, static_cast<PetscInt>(Tag::numVars_))
    {
        static_assert(Tag::numVars_ > 0, "Tag::numVars_ must be positive.");
    }

    /**
     * @brief 兼容旧 façade 的 replicated pure-geometry 语义。
     *
     * 生产 `CpGridCore` 只预计算 owned+ghost；旧 `CpGrid<Tag>` 在 base setup
     * 完成后继续物化完整 Mesh 的 cell/face 纯几何，保持历史远端只读行为。
     */
    void setup()
    {
        CpGridCore::setup();
        for (std::size_t storage = 0; storage < mesh().cellCount(); ++storage)
            cacheCellPureGeometry_(mesh().cellByStorageIndex(storage));
    }
};

} // namespace MPMC
