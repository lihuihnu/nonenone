/**
 * @file structuredgrid.hpp
 * @brief 基于 PETSc DMDA 的规则结构网格及几何/并行访问接口。
 */
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <structuredgrid/gridenums.hpp>
#include <structuredgrid/dmdacontainer.hpp>
#include <structuredgrid/structuredgrid_components.hpp>
#include <structuredgrid/structuredgrid_legacy_geometry.hpp>
#include <structuredgrid/structuredgrid_views.hpp>

namespace MPMC
{

/**
 * @brief 使用物理尺寸定义规则笛卡尔网格范围。
 */
struct GridExtent
{
    double Lx_{0.0};
    double Ly_{0.0};
    double Lz_{0.0};
};

/**
 * @brief 使用均匀单元尺寸定义规则笛卡尔网格。
 */
struct CellSize
{
    double dx_{0.0};
    double dy_{0.0};
    double dz_{0.0};
};

/**
 * @brief 结构化网格内部面的邻接信息。
 *
 * 结构化内部面由邻接单元全局索引和相对于当前单元的六面方向唯一描述。
 */
struct StructuredFace
{
    int neighborCell{-1};
    StructuredFaceOrder direction{StructuredFaceOrder::Count};

    [[nodiscard]] int neighbor() const noexcept
    {
        return neighborCell;
    }
};

/**
 * @brief 最多容纳六个结构网格面的无堆分配小型容器。
 *
 * 该类型提供与只读 `std::vector` 相同的常用遍历接口，但存储直接位于对象内，
 * 避免残差/Jacobian 装配时为每个单元反复申请和释放内存。
 */
class StructuredFaceList final
{
  public:
    using Storage = std::array<StructuredFace, 6>;
    using const_iterator = Storage::const_iterator;
    using size_type = Storage::size_type;

    [[nodiscard]] const_iterator begin() const noexcept
    {
        return storage_.begin();
    }

    [[nodiscard]] const_iterator end() const noexcept
    {
        return storage_.begin() + static_cast<std::ptrdiff_t>(size_);
    }

    [[nodiscard]] size_type size() const noexcept
    {
        return size_;
    }

    [[nodiscard]] bool empty() const noexcept
    {
        return size_ == 0;
    }

    void push_back(StructuredFace face)
    {
        if (size_ >= storage_.size())
            throw std::logic_error("StructuredFaceList capacity exceeded.");
        storage_[size_++] = face;
    }

  private:
    Storage storage_{};
    size_type size_{0};
};

/**
 * @brief 基于 PETSc DMDA 的三维结构化网格封装。
 *
 * StructuredGrid 负责：
 *
 * - 保存规则/分层笛卡尔网格几何信息；
 * - 注册网格自身几何/岩石布局，并允许上层显式注册模型自由度布局；
 * - 创建 PETSc 全局向量、局部向量和矩阵；
 * - 由 `dx/dy/dz` 按需计算体积、面面积和半程距离等可推导几何量；
 * - 仅为 active 与岩石属性保留 PETSc Vec，并提供 DMDA 数组访问接口；
 * - 提供相邻单元、传递系数、扩散系数和重力高度差等几何计算。
 *
 * Grid 核心不理解压力、相、组分、EOS 或 AD。主未知量、phase-state 与
 * 其它模型布局由上层 backend 通过 registerLayout() 显式提供。
 *
 * @note 本类直接管理多个 PETSc Vec 资源，因此禁止拷贝和移动。
 */
class StructuredGridCore
    : private detail::StructuredGridGeometryState,
      private detail::StructuredGridPetscLayoutState,
      private detail::StructuredGridActivityStorage,
      private detail::StructuredGridRockStorage
{
  public:
    using Index = int;
    using Coordinate = std::array<double, 3>;
    using RockWriteView = StructuredGridRockWriteView<StructuredGridCore>;
    using RockReadView = StructuredGridRockReadView<StructuredGridCore>;
    using AssemblyAccess = StructuredGridAssemblyAccess<StructuredGridCore>;

    static_assert(std::is_same_v<PetscScalar, double>,
                  "StructuredGrid requires PETSc configured with real double PetscScalar.");


    // G2 兼容层：真实状态已移动到职责化内部组件。旧字段名暂时通过
    // using 继续暴露，保证现有外部源码无需同步迁移；项目内部应使用 G1 正式 API。
    using detail::StructuredGridGeometryState::Lx_;
    using detail::StructuredGridGeometryState::Ly_;
    using detail::StructuredGridGeometryState::Lz_;
    using detail::StructuredGridGeometryState::nx_;
    using detail::StructuredGridGeometryState::ny_;
    using detail::StructuredGridGeometryState::nz_;
    using detail::StructuredGridGeometryState::dx;
    using detail::StructuredGridGeometryState::dy;
    using detail::StructuredGridGeometryState::dz;
    using detail::StructuredGridGeometryState::activeVec;

    using detail::StructuredGridGeometryState::physicalExtent;
    using detail::StructuredGridGeometryState::dimensions;
    using detail::StructuredGridGeometryState::globalToIJK;
    using detail::StructuredGridGeometryState::ijkToGlobal;
    using detail::StructuredGridGeometryState::contains;
    using detail::StructuredGridGeometryState::cellSize;

    using detail::StructuredGridPetscLayoutState::xl_;
    using detail::StructuredGridPetscLayoutState::yl_;
    using detail::StructuredGridPetscLayoutState::zl_;
    using detail::StructuredGridPetscLayoutState::nxl_;
    using detail::StructuredGridPetscLayoutState::nyl_;
    using detail::StructuredGridPetscLayoutState::nzl_;
    using detail::StructuredGridPetscLayoutState::gxl_;
    using detail::StructuredGridPetscLayoutState::gyl_;
    using detail::StructuredGridPetscLayoutState::gzl_;
    using detail::StructuredGridPetscLayoutState::gnxl_;
    using detail::StructuredGridPetscLayoutState::gnyl_;
    using detail::StructuredGridPetscLayoutState::gnzl_;

    using detail::StructuredGridPetscLayoutState::registerLayout;
    using detail::StructuredGridPetscLayoutState::createGlobalVector;
    using detail::StructuredGridPetscLayoutState::borrowLocalVector;
    using detail::StructuredGridPetscLayoutState::restoreLocalVector;
    using detail::StructuredGridPetscLayoutState::globalToLocal;
    using detail::StructuredGridPetscLayoutState::vecGetArray;
    using detail::StructuredGridPetscLayoutState::vecRestoreArray;
    using detail::StructuredGridPetscLayoutState::vecGetArrayRead;
    using detail::StructuredGridPetscLayoutState::vecRestoreArrayRead;
    using detail::StructuredGridPetscLayoutState::cellIndices;
    using detail::StructuredGridPetscLayoutState::ownedRegion;
    using detail::StructuredGridPetscLayoutState::isSetup;
    using detail::StructuredGridPetscLayoutState::dm;
    using detail::StructuredGridPetscLayoutState::createMatrix;
    using detail::StructuredGridPetscLayoutState::communicator;

    using detail::StructuredGridActivityStorage::active;

    using detail::StructuredGridRockStorage::porosityVector;
    using detail::StructuredGridRockStorage::permeabilityVector;
    using detail::StructuredGridRockStorage::localPermeabilityVector;
    using detail::StructuredGridRockStorage::localPorosityVector;
    using detail::StructuredGridRockStorage::permArr;
    using detail::StructuredGridRockStorage::refPoroArr;

  public:
    /**
     * @brief 使用统一单元尺寸或总体物理范围构造网格。
     *
     * @tparam T 必须为 CellSize 或 GridExtent。
     * @param[in] nx x 方向单元数量。
     * @param[in] ny y 方向单元数量。
     * @param[in] nz z 方向单元数量。
     * @param[in] param 单元尺寸或网格总体范围。
     * @param[in] comm 网格所属 MPI communicator；调用方须保证其生命周期覆盖网格。
     *
     * @throws std::invalid_argument 当网格尺寸或物理长度非法时抛出。
     */
    template <typename T>
    StructuredGridCore(int nx,
                   int ny,
                   int nz,
                   const T &param,
                   MPI_Comm comm = PETSC_COMM_WORLD)
        : detail::StructuredGridGeometryState(nx, ny, nz),
          detail::StructuredGridPetscLayoutState(nx, ny, nz, comm)
    {
        static_assert(std::is_same_v<T, CellSize> ||
                          std::is_same_v<T, GridExtent>,
                      "StructuredGrid 参数必须为 CellSize 或 GridExtent");

        validateGridCounts_();
        if constexpr (std::is_same_v<T, CellSize>)
        {
            validatePositiveLength_(param.dx_, "dx");
            validatePositiveLength_(param.dy_, "dy");
            validatePositiveLength_(param.dz_, "dz");

            dx.assign(static_cast<std::size_t>(nx_), param.dx_);
            dy.assign(static_cast<std::size_t>(ny_), param.dy_);
            dz.assign(static_cast<std::size_t>(nz_), param.dz_);

            Lx_ = param.dx_ * nx_;
            Ly_ = param.dy_ * ny_;
            Lz_ = param.dz_ * nz_;
        }
        else
        {
            validatePositiveLength_(param.Lx_, "Lx");
            validatePositiveLength_(param.Ly_, "Ly");
            validatePositiveLength_(param.Lz_, "Lz");

            Lx_ = param.Lx_;
            Ly_ = param.Ly_;
            Lz_ = param.Lz_;

            dx.assign(static_cast<std::size_t>(nx_), Lx_ / nx_);
            dy.assign(static_cast<std::size_t>(ny_), Ly_ / ny_);
            dz.assign(static_cast<std::size_t>(nz_), Lz_ / nz_);
        }
    }

    /**
     * @brief 使用三个方向独立的均匀尺寸或非均匀尺寸数组构造网格。
     *
     * 每个方向参数可以是 double，也可以是 std::vector<double>。
     * double 表示该方向所有单元采用相同尺寸；vector 表示逐单元尺寸。
     *
     * @tparam T1 x 方向尺寸描述类型。
     * @tparam T2 y 方向尺寸描述类型。
     * @tparam T3 z 方向尺寸描述类型。
     * @param[in] nx x 方向单元数量。
     * @param[in] ny y 方向单元数量。
     * @param[in] nz z 方向单元数量。
     * @param[in] dxInput x 方向单元尺寸。
     * @param[in] dyInput y 方向单元尺寸。
     * @param[in] dzInput z 方向单元尺寸。
     * @param[in] comm 网格所属 MPI communicator；调用方须保证其生命周期覆盖网格。
     *
     * @throws std::invalid_argument 当尺寸数组长度不匹配或存在非正长度时抛出。
     */
    template <typename T1, typename T2, typename T3>
    StructuredGridCore(int nx,
                   int ny,
                   int nz,
                   const T1 &dxInput,
                   const T2 &dyInput,
                   const T3 &dzInput,
                   MPI_Comm comm = PETSC_COMM_WORLD)
        : detail::StructuredGridGeometryState(nx, ny, nz),
          detail::StructuredGridPetscLayoutState(nx, ny, nz, comm)
    {
        validateGridCounts_();
        Lx_ = initializeAxis_(dxInput, nx_, dx, "x");
        Ly_ = initializeAxis_(dyInput, ny_, dy, "y");
        Lz_ = initializeAxis_(dzInput, nz_, dz, "z");
    }

    StructuredGridCore(const StructuredGridCore &) = delete;
    StructuredGridCore &operator=(const StructuredGridCore &) = delete;
    StructuredGridCore(StructuredGridCore &&) = delete;
    StructuredGridCore &operator=(StructuredGridCore &&) = delete;

    /**
     * @brief 释放网格持有的 PETSc Vec 资源。
     *
     * 若调用方遗漏了 unmapLocalArrays()，析构阶段会先尝试
     * 结束尚未完成的数组访问，再释放局部几何向量和岩石属性全局向量。
     */
    ~StructuredGridCore() noexcept
    {
        releaseResources_();
    }

    /**
     * @brief 初始化 canonical DMDA、active ghost 状态和可推导几何缓存。
     *
     * 主要步骤：
     *
     * 1. 注册 Grid 自身需要的 1/3 DOF 布局；
     * 2. 获取当前 MPI 进程的 owned/ghost 范围；
     * 3. 按既有累积顺序缓存 O(nz) 垂向中心，用于重力势差；
     * 4. 仅同步 active 标记到长期局部 ghost Vec；
     * 5. 构造当前进程拥有的活动单元全局索引列表。
     *
     * 半程距离、面面积和 bulk volume 不再存为 O(Ncell) PETSc Vec，
     * 而是在使用点由轴宽直接计算。
     *
     * @throws std::logic_error 当重复调用 setup() 时抛出。
     * @throws std::invalid_argument 当 activeVec 尺寸与网格不匹配时抛出。
     *
     * @note 若 activeVec 为空，默认所有单元均为活动单元。
     */
    void setup()
    {
        if (setupComplete_)
        {
            throw std::logic_error("StructuredGrid::setup() 不能重复调用");
        }

        validateGeometry_();
        initializeActivity_();
        cacheVerticalCenters_();
        registerDofLayouts_();

        const DMDARegion owned = dmCont_.layout(1).ownedRegion();
        xl_ = owned.xStart;
        yl_ = owned.yStart;
        zl_ = owned.zStart;
        nxl_ = owned.xCount;
        nyl_ = owned.yCount;
        nzl_ = owned.zCount;

        const DMDARegion ghost = dmCont_.layout(1).ghostRegion();
        gxl_ = ghost.xStart;
        gyl_ = ghost.yStart;
        gzl_ = ghost.zStart;
        gnxl_ = ghost.xCount;
        gnyl_ = ghost.yCount;
        gnzl_ = ghost.zCount;

        Vec activeGlobal = createGlobalVector(1);
        auto activeArray = vecGetArray<1>(activeGlobal);

        cellIndices_.clear();
        cellIndices_.reserve(static_cast<std::size_t>(nxl_) *
                             static_cast<std::size_t>(nyl_) *
                             static_cast<std::size_t>(nzl_));

        for (int k = zl_; k < zl_ + nzl_; ++k)
        {
            for (int j = yl_; j < yl_ + nyl_; ++j)
            {
                for (int i = xl_; i < xl_ + nxl_; ++i)
                {
                    const int globalIndex =
                        index_.ijkToGlobalUnchecked(i, j, k);
                    activeArray[k][j][i][0] = activeVec[globalIndex];
                    if (activeArray[k][j][i][0] != 0.0)
                        cellIndices_.push_back(globalIndex);
                }
            }
        }

        vecRestoreArray<1>(activeGlobal, activeArray);
        active = createPersistentLocal_(1, activeGlobal);
        cacheLocalActivity_();
        destroyVec_(activeGlobal);

        // activeVec 只用于 setup() 阶段，可释放其全局复制以降低 MPI 大网格内存。
        activeVec.clear();
        activeVec.shrink_to_fit();

        // dx/dy/dz 在 setup() 后继续保留，供 cellSize() 和 Peaceman 井指数计算读取。
        setupComplete_ = true;
    }

    /**
     * @brief 获取带 ghost 的局部孔隙度/渗透率数组。
     *
     * localPermeabilityVector 和 localPorosityVector 通过 DMGetLocalVector() 从 PETSc
     * 局部向量池临时借用，并在 unmapLocalArrays() 中归还。
     *
     * @warning 必须与 unmapLocalArrays() 成对调用。
     */
    void mapLocalArrays()
    {
        ensureRockPropertiesInitialized_();
        beginArrayAccess_();

        localPermeabilityVector = borrowLocalVector(3, permeabilityVector);
        localPorosityVector = borrowLocalVector(1, porosityVector);

        permArr = vecGetArray<3>(localPermeabilityVector);
        refPoroArr = vecGetArray<1>(localPorosityVector);
    }

    /**
     * @brief 恢复局部数组映射，并将临时局部岩石属性向量归还 PETSc。
     *
     * 必须恢复 mapLocalArrays() 中取得的同一 Vec/数组指针对，并归还临时 local Vec。
     */
    void unmapLocalArrays()
    {
        if (!arraysMapped_)
            throw std::logic_error(
                "StructuredGrid::unmapLocalArrays() 没有对应的 mapLocalArrays()");

        vecRestoreArray<3>(localPermeabilityVector, permArr);
        vecRestoreArray<1>(localPorosityVector, refPoroArr);

        clearRockArrayPointers_();

        dmCont_.layout(3).restoreLocalVector(localPermeabilityVector);
        dmCont_.layout(1).restoreLocalVector(localPorosityVector);

        arraysMapped_ = false;
    }

    /**
     * @brief 设置 setup() 阶段使用的全局活动单元输入。
     *
     * 输入值的有限性、尺寸和 0/1 归一化仍由 setup() 中的既有校验完成，
     * 因此本接口不改变原 activeVec 的验证时机。
     */
    void setActivityInput(std::vector<double> activity)
    {
        if (setupComplete_)
            throw std::logic_error(
                "StructuredGrid activity input must be set before setup().");
        activeVec = std::move(activity);
    }

    /** @brief 返回拥有型三向渗透率全局 Vec 的只读非拥有句柄副本。 */
    [[nodiscard]] Vec permeabilityVectorHandle() const noexcept
    {
        return permeabilityVector;
    }

    /** @brief 返回拥有型孔隙度全局 Vec 的只读非拥有句柄副本。 */
    [[nodiscard]] Vec porosityVectorHandle() const noexcept
    {
        return porosityVector;
    }

    /** @brief 建立全局岩石属性的可写 RAII 视图。 */
    [[nodiscard]] RockWriteView rockWriteView()
    {
        return RockWriteView(*this);
    }

    /** @brief 建立全局岩石属性的只读 RAII 视图。 */
    [[nodiscard]] RockReadView rockReadView() const
    {
        return RockReadView(*this);
    }

    /**
     * @brief 建立装配期 ghosted 几何/岩石数组访问守卫。
     *
     * 其生命周期内现有 porosity()/cellVolume()/transmissibility()/
     * gravityPotentialDifference() 等接口保持可用；析构时自动恢复所有映射。
     */
    [[nodiscard]] AssemblyAccess assemblyAccess()
    {
        return AssemblyAccess(*this);
    }

    /**
     * @brief 创建并接管标准的三向渗透率和标量孔隙度全局向量。
     *
     * 该接口替代调用方分别创建并直接赋值两个拥有型 Vec 的易错写法。
     * 已初始化时重复调用会抛出异常，避免覆盖句柄并泄漏 PETSc 资源。
     */
    void initializeRockProperties()
    {
        if (!setupComplete_)
            throw std::logic_error(
                "StructuredGrid::initializeRockProperties requires setup() first.");
        if (permeabilityVector != nullptr || porosityVector != nullptr)
            throw std::logic_error(
                "StructuredGrid rock properties are already initialized.");

        permeabilityVector = createGlobalVector(3);
        porosityVector = createGlobalVector(1);
    }

    /** @brief 查询标准岩石属性向量是否均已初始化。 */
    [[nodiscard]] bool hasRockProperties() const noexcept
    {
        return permeabilityVector != nullptr && porosityVector != nullptr;
    }

    /**
     * @brief 判断全局单元是否为当前网格中的活动单元。
     *
     * @note 对当前 rank 拥有单元及其 stencil ghost 邻居有效。
     */
    [[nodiscard]] bool isActive(int globalIndex) const
    {
        if (!setupComplete_)
        {
            throw std::logic_error(
                "StructuredGrid::isActive requires setup() first.");
        }
        return isActiveCell_(globalIndex);
    }

    /**
     * @brief 获取 StructuredFace 对应的邻接单元全局索引。
     */
    [[nodiscard]] int neighbors(const StructuredFace &face) const noexcept
    {
        return face.neighborCell;
    }

    /**
     * @brief 获取指定活动单元的所有活动内部邻接面。
     *
     * 几何越界邻居和 inactive 邻居均不会返回，因此 active/inactive 接口被视为
     * 无流边界。共享内部面仍会从两个活动控制体分别遍历，以保持当前残差与
     * Jacobian 的双侧装配方式。
     */
    [[nodiscard]] StructuredFaceList faces(int globalIndex) const
    {
        if (!setupComplete_)
        {
            throw std::logic_error(
                "StructuredGrid::faces requires setup() first.");
        }

        if (!onProcess(globalIndex))
        {
            throw std::out_of_range(
                "StructuredGrid::faces expects a cell owned by the current MPI rank.");
        }

        StructuredFaceList result;

        if (!isActiveCell_(globalIndex))
        {
            return result;
        }

        const auto [i, j, k] = index_.globalToIJK(globalIndex);

        const auto append = [&](int ni,
                                int nj,
                                int nk,
                                StructuredFaceOrder direction)
        {
            if (!index_.contains(ni, nj, nk))
            {
                return;
            }

            const int neighbor = index_.ijkToGlobalUnchecked(ni, nj, nk);
            if (!isActiveCell_(neighbor))
            {
                return;
            }

            result.push_back({neighbor, direction});
        };

        append(i - 1, j, k, StructuredFaceOrder::XMinus);
        append(i + 1, j, k, StructuredFaceOrder::XPlus);
        append(i, j - 1, k, StructuredFaceOrder::YMinus);
        append(i, j + 1, k, StructuredFaceOrder::YPlus);
        append(i, j, k - 1, StructuredFaceOrder::ZMinus);
        append(i, j, k + 1, StructuredFaceOrder::ZPlus);

        return result;
    }

    /**
     * @brief 计算两个相邻单元中心的重力势差项 `g * Δz`。
     *
     * @param[in] i 第一个单元全局索引。
     * @param[in] j 第二个单元全局索引。
     * @return `9.80665 * (z_i - z_j)`。
     *
     * @pre 必须已建立装配期数组访问；垂向中心来自 O(nz) 一维缓存。
     */
    [[nodiscard]] double gravityPotentialDifference(int i, int j) const
    {
        return gravityPotentialDifference(
            i, StructuredFace{j, directionBetween_(i, j)});
    }

    /**
     * @brief 使用 faces() 已知的方向信息计算重力势差，避免重复解析邻接方向。
     */
    [[nodiscard]] double gravityPotentialDifference(
        int cell,
        const StructuredFace &face) const
    {
        ensureArraysMapped_();
        validateFace_(cell, face);
        if (structuredFaceAxis(face.direction) != 2)
            return 0.0;

        const auto [i, j, k] = index_.globalToIJKUnchecked(cell);
        const int neighborK = k + structuredFaceSign(face.direction);
        return (verticalCenter_(k) - verticalCenter_(neighborK)) *
               gravityAcceleration_;
    }

    /**
     * @brief 计算相邻单元之间的两点通量传递系数。
     *
     * 两侧半传递系数按调和方式串联：
     * `T = 1 / (1/T_i + 1/T_j)`。
     *
     * @pre 必须通过装配期数组访问获取 ghosted 渗透率；几何量按需计算。
     */
    [[nodiscard]] double transmissibility(int i, int j) const
    {
        return transmissibility(i, StructuredFace{j, directionBetween_(i, j)});
    }

    /**
     * @brief 使用 faces() 已知的方向计算 TPFA 传递系数。
     */
    [[nodiscard]] double transmissibility(
        int cell,
        const StructuredFace &face) const
    {
        ensureArraysMapped_();
        validateFace_(cell, face);

        const auto [i, j, k] = index_.globalToIJKUnchecked(cell);
        const auto [ni, nj, nk] = neighborCoordinates_(i, j, k, face.direction);
        const double halfTrans1 = computeHalfTrans_(i, j, k, face.direction);
        const double halfTrans2 = computeHalfTrans_(
            ni, nj, nk, oppositeStructuredFace(face.direction));

        if (halfTrans1 <= transmissibilityEpsilon_ ||
            halfTrans2 <= transmissibilityEpsilon_)
        {
            return 0.0;
        }

        const double result =
            1.0 / (1.0 / halfTrans1 + 1.0 / halfTrans2);
        if (!std::isfinite(result) || result < 0.0)
            throw std::runtime_error(
                "StructuredGrid computed an invalid transmissibility.");
        return result;
    }

    /**
     * @brief 获取指定单元体积。
     */
    [[nodiscard]] double cellVolume(int globalIndex) const
    {
        ensureArraysMapped_();
        const auto [i, j, k] = index_.globalToIJK(globalIndex);
        return detail::structuredCellVolume({
            dx[static_cast<std::size_t>(i)],
            dy[static_cast<std::size_t>(j)],
            dz[static_cast<std::size_t>(k)]});
    }

    /**
     * @brief 获取指定全局单元的参考孔隙度。
     */
    [[nodiscard]] double porosity(int globalIndex) const
    {
        ensureArraysMapped_();
        const auto [i, j, k] = index_.globalToIJK(globalIndex);
        const double value = refPoroArr[k][j][i][0];
        if (!std::isfinite(value) || value < 0.0 || value > 1.0)
            throw std::runtime_error(
                "StructuredGrid porosity must be finite and within [0, 1].");
        return value;
    }

    /**
     * @brief 判断指定全局单元是否由当前 MPI 进程拥有。
     */
    [[nodiscard]] bool onProcess(int globalIndex) const noexcept
    {
        if (!index_.contains(globalIndex))
            return false;
        const auto [i, j, k] = index_.globalToIJKUnchecked(globalIndex);
        return xl_ <= i && i < xl_ + nxl_ && yl_ <= j &&
               j < yl_ + nyl_ && zl_ <= k && k < zl_ + nzl_;
    }

  private:
    static constexpr double gravityAcceleration_{9.80665};
    static constexpr double transmissibilityEpsilon_{1.0e-30};

    /**
     * @brief 注册 StructuredGrid 当前已知的数据布局。
     */
    void registerDofLayouts_()
    {
        // Grid 核心只注册 active/porosity 的 1-DOF 与 permeability 的 3-DOF。
        // 可推导 Cartesian 几何不再需要 6/18-DOF DMDA；模型主未知量、
        // phase-state 和 AD 辅助布局仍由上层模型显式注册。
        dmCont_.registerLayout(1);
        dmCont_.registerLayout(3);
    }

    /** @brief 开始唯一受支持的 ghosted 数组映射周期。 */
    void beginArrayAccess_()
    {
        if (!setupComplete_)
            throw std::logic_error(
                "StructuredGrid 必须先调用 setup() 才能访问网格数组");
        if (arraysMapped_)
            throw std::logic_error(
                "StructuredGrid 已存在尚未恢复的 PETSc 数组访问");
        arraysMapped_ = true;
    }

    /** @brief 防止在 PETSc 数组映射周期之外解引用空数组指针。 */
    void ensureArraysMapped_() const
    {
        if (!arraysMapped_)
            throw std::logic_error(
                "StructuredGrid assembly access requires mapLocalArrays().");
    }

    /**
     * @brief 检查孔隙度和渗透率全局向量是否已经由调用方配置。
     */
    void ensureRockPropertiesInitialized_() const
    {
        if (permeabilityVector == nullptr || porosityVector == nullptr)
        {
            throw std::logic_error(
                "StructuredGrid::permeabilityVector 和 porosityVector 必须先初始化");
        }
    }

    /**
     * @brief 清除岩石属性数组指针。
     */
    void clearRockArrayPointers_() noexcept
    {
        permArr = nullptr;
        refPoroArr = nullptr;
    }

    /**
     * @brief 将 ghosted active Vec 缓存为紧凑的本地活动标记。
     *
     * 这样 faces() 可以在不反复调用 PETSc 数组映射接口的情况下判断邻居是否
     * 活动，同时仅保存当前 rank 的 owned+ghost 区域，而不是复制整个全局 activeVec。
     */
    void cacheLocalActivity_()
    {
        localActiveMask_.assign(
            static_cast<std::size_t>(gnxl_) * static_cast<std::size_t>(gnyl_) *
                static_cast<std::size_t>(gnzl_),
            0);

        auto activeArray = vecGetArray<1>(active);
        for (int k = gzl_; k < gzl_ + gnzl_; ++k)
        {
            for (int j = gyl_; j < gyl_ + gnyl_; ++j)
            {
                for (int i = gxl_; i < gxl_ + gnxl_; ++i)
                {
                    localActiveMask_[localGhostOffset_(i, j, k)] =
                        activeArray[k][j][i][0] != 0.0 ? 1 : 0;
                }
            }
        }
        vecRestoreArray<1>(active, activeArray);
    }

    [[nodiscard]] std::size_t
    localGhostOffset_(int i, int j, int k) const
    {
        if (i < gxl_ || i >= gxl_ + gnxl_ || j < gyl_ ||
            j >= gyl_ + gnyl_ || k < gzl_ || k >= gzl_ + gnzl_)
        {
            throw std::out_of_range(
                "StructuredGrid ghost lookup is outside the local ghost region.");
        }

        const std::size_t ii = static_cast<std::size_t>(i - gxl_);
        const std::size_t jj = static_cast<std::size_t>(j - gyl_);
        const std::size_t kk = static_cast<std::size_t>(k - gzl_);
        return (kk * static_cast<std::size_t>(gnyl_) + jj) *
                   static_cast<std::size_t>(gnxl_) +
               ii;
    }

    [[nodiscard]] bool isActiveCell_(int globalIndex) const
    {
        const auto [i, j, k] = index_.globalToIJK(globalIndex);
        return localActiveMask_.at(localGhostOffset_(i, j, k)) != 0;
    }

    /**
     * @brief 计算正交结构网格当前单元一侧的 TPFA 半传递系数。
     *
     * 对从单元中心指向面中心的向量 `d`、面法向面积向量 `A_n` 与对应
     * 主方向渗透率 `K`，使用 `T_half = K |A_n . d| / |d|^2`。相邻两侧
     * 半传递系数在 `transmissibility()` 中按调和平均串联。
     */
    [[nodiscard]] static std::array<int, 3> neighborCoordinates_(
        int i,
        int j,
        int k,
        StructuredFaceOrder direction)
    {
        switch (direction)
        {
        case StructuredFaceOrder::XMinus: --i; break;
        case StructuredFaceOrder::XPlus:  ++i; break;
        case StructuredFaceOrder::YMinus: --j; break;
        case StructuredFaceOrder::YPlus:  ++j; break;
        case StructuredFaceOrder::ZMinus: --k; break;
        case StructuredFaceOrder::ZPlus:  ++k; break;
        case StructuredFaceOrder::Count:
            throw std::invalid_argument("Invalid StructuredFace direction.");
        }
        return {i, j, k};
    }

    /** @brief 由两个相邻全局索引确定从第一个单元指向第二个单元的方向。 */
    [[nodiscard]] StructuredFaceOrder directionBetween_(int cell, int neighbor) const
    {
        const auto [i, j, k] = index_.globalToIJK(cell);
        const auto [ni, nj, nk] = index_.globalToIJK(neighbor);

        if (ni == i - 1 && nj == j && nk == k) return StructuredFaceOrder::XMinus;
        if (ni == i + 1 && nj == j && nk == k) return StructuredFaceOrder::XPlus;
        if (ni == i && nj == j - 1 && nk == k) return StructuredFaceOrder::YMinus;
        if (ni == i && nj == j + 1 && nk == k) return StructuredFaceOrder::YPlus;
        if (ni == i && nj == j && nk == k - 1) return StructuredFaceOrder::ZMinus;
        if (ni == i && nj == j && nk == k + 1) return StructuredFaceOrder::ZPlus;

        throw std::invalid_argument(
            "StructuredGrid transmissibility requires face-adjacent cells.");
    }

    /** @brief 校验 StructuredFace 的方向与邻居索引彼此一致。 */
    void validateFace_(int cell, const StructuredFace &face) const
    {
        if (!index_.contains(cell) || !index_.contains(face.neighborCell) ||
            !isValidStructuredFace(face.direction))
            throw std::out_of_range("StructuredGrid face is invalid.");

        const auto [i, j, k] = index_.globalToIJKUnchecked(cell);
        const auto neighbor = neighborCoordinates_(i, j, k, face.direction);
        if (!index_.contains(neighbor[0], neighbor[1], neighbor[2]) ||
            index_.ijkToGlobalUnchecked(neighbor[0], neighbor[1], neighbor[2]) !=
                face.neighborCell)
            throw std::invalid_argument(
                "StructuredGrid face direction does not match its neighbor.");
    }

    [[nodiscard]] double computeHalfTrans_(
        int i,
        int j,
        int k,
        StructuredFaceOrder direction) const
    {
        const int axis = structuredFaceAxis(direction);
        if (axis < 0)
            throw std::invalid_argument("Invalid StructuredFace direction.");

        const double permeabilityValue = permArr[k][j][i][axis];
        if (!std::isfinite(permeabilityValue) || permeabilityValue < 0.0)
            throw std::runtime_error(
                "StructuredGrid permeability must be finite and non-negative.");

        const std::array<double, 3> widths{
            dx[static_cast<std::size_t>(i)],
            dy[static_cast<std::size_t>(j)],
            dz[static_cast<std::size_t>(k)]};
        const double distanceValue =
            detail::structuredFaceHalfDistance(widths, direction);
        const double normalComponent =
            detail::structuredFaceNormalComponent(widths, direction);

        if (!std::isfinite(normalComponent) ||
            !std::isfinite(distanceValue) || distanceValue <= 0.0)
            throw std::runtime_error(
                "StructuredGrid encountered invalid face geometry.");

        return permeabilityValue * std::abs(normalComponent) / distanceValue;
    }

    /**
     * @brief 销毁一个由本类拥有的 PETSc Vec。
     */
    static void destroyVec_(Vec &vec) noexcept
    {
        if (vec == nullptr)
        {
            return;
        }
        PetscCallAbort(PETSC_COMM_SELF, VecDestroy(&vec));
    }

    /** @brief 析构时回收尚未结束的 ghosted 数组映射。 */
    void restoreOutstandingArrays_() noexcept
    {
        if (!arraysMapped_)
            return;

        if (permArr != nullptr && localPermeabilityVector != nullptr)
            vecRestoreArray<3>(localPermeabilityVector, permArr);
        if (refPoroArr != nullptr && localPorosityVector != nullptr)
            vecRestoreArray<1>(localPorosityVector, refPoroArr);

        if (localPermeabilityVector != nullptr && dmCont_.contains(3))
            dmCont_.layout(3).restoreLocalVector(localPermeabilityVector);
        if (localPorosityVector != nullptr && dmCont_.contains(1))
            dmCont_.layout(1).restoreLocalVector(localPorosityVector);

        permArr = nullptr;
        refPoroArr = nullptr;
        arraysMapped_ = false;
    }

    /**
     * @brief 释放 StructuredGrid 持有的 PETSc 资源。
     */
    void releaseResources_() noexcept
    {
        restoreOutstandingArrays_();

        destroyVec_(active);

        destroyVec_(porosityVector);
        destroyVec_(permeabilityVector);

        // 正常情况下临时局部向量已经在 unmapLocalArrays() 中归还。
        if (localPermeabilityVector != nullptr && dmCont_.contains(3))
        {
            dmCont_.layout(3).restoreLocalVector(localPermeabilityVector);
        }
        if (localPorosityVector != nullptr && dmCont_.contains(1))
        {
            dmCont_.layout(1).restoreLocalVector(localPorosityVector);
        }
    }
 };

/**
 * @brief 兼容旧 `StructuredGrid<Tag>` 源码接口的薄 façade。
 *
 * Grid 核心不再读取模型 Tag。该兼容层仅在 `setup()` 前按旧顺序注册
 * Tag 提供的布局，使历史调用保持原有 canonical DMDA 和布局注册行为。
 * 新代码应优先直接使用 `StructuredGridCore`，由模型/backend 显式注册布局。
 */
template <class Tag>
class StructuredGrid final
    : public StructuredGridCore,
      private detail::StructuredGridLegacyGeometryCompatibility
{
  public:
    using detail::StructuredGridLegacyGeometryCompatibility::distance;
    using detail::StructuredGridLegacyGeometryCompatibility::areaNormal;
    using detail::StructuredGridLegacyGeometryCompatibility::center;
    using detail::StructuredGridLegacyGeometryCompatibility::volume;
    using detail::StructuredGridLegacyGeometryCompatibility::distanceArr;
    using detail::StructuredGridLegacyGeometryCompatibility::areaNormalArr;
    using detail::StructuredGridLegacyGeometryCompatibility::centerArr;
    using detail::StructuredGridLegacyGeometryCompatibility::volumeArr;

    template <typename T>
    StructuredGrid(int nx,
                   int ny,
                   int nz,
                   const T &param,
                   MPI_Comm comm = PETSC_COMM_WORLD)
        : StructuredGridCore(nx, ny, nz, param, comm)
    {
    }

    template <typename T1, typename T2, typename T3>
    StructuredGrid(int nx,
                   int ny,
                   int nz,
                   const T1 &dxInput,
                   const T2 &dyInput,
                   const T3 &dzInput,
                   MPI_Comm comm = PETSC_COMM_WORLD)
        : StructuredGridCore(nx, ny, nz, dxInput, dyInput, dzInput, comm)
    {
    }

    ~StructuredGrid() noexcept
    {
        restoreLegacyGeometryArrays_(*this);
        destroyLegacyGeometry_();
    }

    void setup()
    {
        static_assert(Tag::numVars_ > 0, "Tag::numVars_ must be positive.");
        static_assert(Tag::numPhaseState_ > 0, "Tag::numPhaseState_ must be positive.");

        // 完整保留旧 façade 的 layout 注册顺序。生产 Core 自身只需要 1/3 DOF。
        registerLayout(Tag::numVars_);
        if constexpr (!std::is_same_v<typename Tag::ValueType, double>)
            registerLayout(1 + Tag::numVars_);
        registerLayout(Tag::numPhaseState_);
        registerLayout(1);
        registerLayout(3);
        registerLayout(6);
        registerLayout(18);

        StructuredGridCore::setup();
        initializeLegacyGeometry_(*this);
    }

    /**
     * @brief 旧 façade 手工数组接口：在 Core 岩石映射之外继续映射历史几何 Vec。
     */
    void mapLocalArrays()
    {
        StructuredGridCore::mapLocalArrays();
        mapLegacyGeometryArrays_(*this);
    }

    /** @brief 恢复旧 façade 的几何数组，然后恢复 Core 岩石数组。 */
    void unmapLocalArrays()
    {
        restoreLegacyGeometryArrays_(*this);
        StructuredGridCore::unmapLocalArrays();
    }


};

} // namespace MPMC
