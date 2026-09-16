/**
 * @file dmdaelem.cpp
 * @brief StructuredGrid 单元、面和邻接关系辅助的实现。
 */
#include <structuredgrid/dmdaelem.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <utility>

namespace MPMC
{
namespace
{

void validateConfiguration(int nx,
                           int ny,
                           int nz,
                           int dof,
                           int stencilWidth,
                           MPI_Comm comm)
{
    if (nx <= 0 || ny <= 0 || nz <= 0)
    {
        throw std::invalid_argument(
            "DMDAElem requires positive grid dimensions.");
    }
    if (dof <= 0)
    {
        throw std::invalid_argument("DMDAElem requires a positive dof.");
    }
    if (stencilWidth < 0)
    {
        throw std::invalid_argument(
            "DMDAElem requires a non-negative stencil width.");
    }
    if (comm == MPI_COMM_NULL)
    {
        throw std::invalid_argument(
            "DMDAElem requires a valid MPI communicator.");
    }
}

void validateDof(int dof)
{
    if (dof <= 0)
    {
        throw std::invalid_argument("DMDAElem requires a positive dof.");
    }
}

/** @brief 为薄层规则网格选择不超过各方向单元数的 MPI 进程拓扑。 */
[[nodiscard]] std::array<int, 3> chooseProcessorGrid(
    int nx,
    int ny,
    int nz,
    MPI_Comm comm)
{
    int processCount = 0;
    if (MPI_Comm_size(comm, &processCount) != MPI_SUCCESS)
        throw std::runtime_error("Unable to query the MPI communicator size.");

    std::array<int, 3> best{0, 0, 0};
    double bestSurface = std::numeric_limits<double>::infinity();
    double bestLongestSide = std::numeric_limits<double>::infinity();

    for (int px = 1; px <= std::min(nx, processCount); ++px)
    {
        if (processCount % px != 0)
            continue;
        const int yzProcesses = processCount / px;
        for (int py = 1; py <= std::min(ny, yzProcesses); ++py)
        {
            if (yzProcesses % py != 0)
                continue;
            const int pz = yzProcesses / py;
            if (pz > nz)
                continue;

            const double lx = static_cast<double>(nx) / px;
            const double ly = static_cast<double>(ny) / py;
            const double lz = static_cast<double>(nz) / pz;
            const double surface = lx * ly + lx * lz + ly * lz;
            const double longestSide = std::max({lx, ly, lz});
            if (surface < bestSurface ||
                (surface == bestSurface && longestSide < bestLongestSide))
            {
                best = {px, py, pz};
                bestSurface = surface;
                bestLongestSide = longestSide;
            }
        }
    }

    if (best[0] == 0)
    {
        throw std::invalid_argument(
            "StructuredGrid requires no more MPI ranks than grid cells.");
    }
    return best;
}

[[nodiscard]] int toInt(PetscInt value)
{
    if (value < static_cast<PetscInt>(std::numeric_limits<int>::min()) ||
        value > static_cast<PetscInt>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error(
            "PETSc integer exceeds the range supported by MPMC grid indices.");
    }
    return static_cast<int>(value);
}

[[nodiscard]] MPI_Comm communicator(DM dm) noexcept
{
    return PetscObjectComm(reinterpret_cast<PetscObject>(dm));
}

void requireDM(DM dm)
{
    if (dm == nullptr)
    {
        throw std::logic_error("DMDAElem does not own a valid PETSc DM.");
    }
}

void requireVec(Vec vec, const char *name)
{
    if (vec == nullptr)
    {
        throw std::invalid_argument(name);
    }
}

} // namespace

DMDAElem::DMDAElem(int nx,
                   int ny,
                   int nz,
                   int dof,
                   int stencilWidth,
                   DMBoundaryType boundaryType,
                   DMDAStencilType stencilType,
                   MPI_Comm comm)
{
    validateConfiguration(nx, ny, nz, dof, stencilWidth, comm);
    const auto processors = chooseProcessorGrid(nx, ny, nz, comm);

    PetscCallAbort(
        comm,
        DMDACreate3d(comm,
                     boundaryType,
                     boundaryType,
                     boundaryType,
                     stencilType,
                     nx,
                     ny,
                     nz,
                     processors[0],
                     processors[1],
                     processors[2],
                     dof,
                     stencilWidth,
                     nullptr,
                     nullptr,
                     nullptr,
                     &dm_));

    PetscCallAbort(comm, DMSetFromOptions(dm_));
    PetscCallAbort(comm, DMSetUp(dm_));

    try
    {
        refreshMetadata_();
        if (dof_ != dof)
        {
            throw std::runtime_error(
                "PETSc changed the requested DMDA dof unexpectedly.");
        }
    }
    catch (...)
    {
        PetscCallAbort(comm, DMDestroy(&dm_));
        throw;
    }
}

DMDAElem::DMDAElem(const DMDAElem &layoutSource, int dof)
{
    requireDM(layoutSource.dm_);
    validateDof(dof);

    const MPI_Comm comm = communicator(layoutSource.dm_);
    PetscCallAbort(
        comm,
        DMDACreateCompatibleDMDA(layoutSource.dm_, dof, &dm_));

    try
    {
        refreshMetadata_();
        if (dof_ != dof)
        {
            throw std::runtime_error(
                "Compatible DMDA did not preserve the requested dof.");
        }
    }
    catch (...)
    {
        PetscCallAbort(comm, DMDestroy(&dm_));
        throw;
    }
}

DMDAElem::DMDAElem(DMDAElem &&other) noexcept
    : dm_(std::exchange(other.dm_, nullptr)),
      dof_(other.dof_),
      dimensions_(other.dimensions_),
      ownedRegion_(other.ownedRegion_),
      ghostRegion_(other.ghostRegion_)
{
    other.reset_();
}

DMDAElem &DMDAElem::operator=(DMDAElem &&other) noexcept
{
    if (this == &other)
    {
        return *this;
    }

    release_();
    dm_ = std::exchange(other.dm_, nullptr);
    dof_ = other.dof_;
    dimensions_ = other.dimensions_;
    ownedRegion_ = other.ownedRegion_;
    ghostRegion_ = other.ghostRegion_;
    other.reset_();
    return *this;
}

DMDAElem::~DMDAElem() noexcept
{
    release_();
}

Vec DMDAElem::createGlobalVector() const
{
    requireDM(dm_);
    Vec vector = nullptr;
    PetscCallAbort(communicator(dm_), DMCreateGlobalVector(dm_, &vector));
    return vector;
}

Vec DMDAElem::borrowLocalVector() const
{
    requireDM(dm_);
    Vec vector = nullptr;
    PetscCallAbort(communicator(dm_), DMGetLocalVector(dm_, &vector));
    return vector;
}

Vec DMDAElem::borrowLocalVector(Vec global) const
{
    requireVec(global,
               "DMDAElem::borrowLocalVector requires a valid global Vec.");
    Vec local = borrowLocalVector();
    globalToLocal(global, local);
    return local;
}

void DMDAElem::restoreLocalVector(Vec &local) const
{
    requireDM(dm_);
    if (local == nullptr)
    {
        return;
    }
    PetscCallAbort(communicator(dm_), DMRestoreLocalVector(dm_, &local));
}

void DMDAElem::globalToLocal(Vec global, Vec local) const
{
    requireDM(dm_);
    requireVec(global, "DMDAElem::globalToLocal requires a valid global Vec.");
    requireVec(local, "DMDAElem::globalToLocal requires a valid local Vec.");

    const MPI_Comm comm = communicator(dm_);
    PetscCallAbort(comm,
                   DMGlobalToLocalBegin(dm_, global, INSERT_VALUES, local));
    PetscCallAbort(comm,
                   DMGlobalToLocalEnd(dm_, global, INSERT_VALUES, local));
}

Mat DMDAElem::createMatrix() const
{
    requireDM(dm_);

    Mat matrix = nullptr;
    const MPI_Comm comm = communicator(dm_);
    PetscCallAbort(comm, DMCreateMatrix(dm_, &matrix));

    // 保留现有 Jacobian 装配策略：允许出现预分配模式之外的新非零位置。
    PetscCallAbort(
        comm,
        MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE));
    return matrix;
}

void DMDAElem::refreshMetadata_()
{
    requireDM(dm_);
    const MPI_Comm comm = communicator(dm_);

    PetscInt xs = 0;
    PetscInt ys = 0;
    PetscInt zs = 0;
    PetscInt xm = 0;
    PetscInt ym = 0;
    PetscInt zm = 0;

    PetscInt gxs = 0;
    PetscInt gys = 0;
    PetscInt gzs = 0;
    PetscInt gxm = 0;
    PetscInt gym = 0;
    PetscInt gzm = 0;

    PetscInt nx = 0;
    PetscInt ny = 0;
    PetscInt nz = 0;
    PetscInt dof = 0;

    PetscCallAbort(
        comm,
        DMDAGetCorners(dm_, &xs, &ys, &zs, &xm, &ym, &zm));
    PetscCallAbort(
        comm,
        DMDAGetGhostCorners(dm_, &gxs, &gys, &gzs, &gxm, &gym, &gzm));
    PetscCallAbort(
        comm,
        DMDAGetInfo(dm_,
                    nullptr,
                    &nx,
                    &ny,
                    &nz,
                    nullptr,
                    nullptr,
                    nullptr,
                    &dof,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr));

    dimensions_ = {toInt(nx), toInt(ny), toInt(nz)};
    dof_ = toInt(dof);
    ownedRegion_ = {toInt(xs), toInt(ys), toInt(zs),
                    toInt(xm), toInt(ym), toInt(zm)};
    ghostRegion_ = {toInt(gxs), toInt(gys), toInt(gzs),
                    toInt(gxm), toInt(gym), toInt(gzm)};
}

void DMDAElem::release_() noexcept
{
    if (dm_ == nullptr)
    {
        return;
    }
    const MPI_Comm comm = communicator(dm_);
    PetscCallAbort(comm, DMDestroy(&dm_));
    reset_();
}

void DMDAElem::reset_() noexcept
{
    dof_ = 0;
    dimensions_ = {};
    ownedRegion_ = {};
    ghostRegion_ = {};
}

} // namespace MPMC
