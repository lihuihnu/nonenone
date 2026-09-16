/**
 * @file dof_layout.cpp
 * @brief CpGrid 单元自由度布局及 PETSc owned/ghost 映射的实现。
 */
#include <cpgrid/dof_layout.hpp>

#include <petscdmshell.h>

#include <numeric>
#include <stdexcept>
#include <vector>

namespace MPMC
{

DofLayout::DofLayout(
    Mesh &mesh,
    PetscInt dofPerCell)
    : mesh_(mesh),
      dofMap_(mesh, dofPerCell)
{
    // mapping 与首个 global->local scatter 使用完全相同的标量索引。
    // 构造阶段只展开一次，避免高 DOF layout 重复生成 O(Nlocal*dof) 临时数组。
    const auto globalIndexValues =
        dofMap_.localToGlobalIndices();
    createMapping_(globalIndexValues);

    try
    {
        createShell_();
        createGlobalToLocalScatter_(globalIndexValues);
    }
    catch (...)
    {
        release_();
        throw;
    }
}

DofLayout::~DofLayout() noexcept
{
    release_();
}

PetscErrorCode DofLayout::context_(
    DM dm,
    DofLayout **layout)
{
    PetscFunctionBegin;

    PetscCheck(
        layout != nullptr,
        PetscObjectComm(
            reinterpret_cast<PetscObject>(dm)),
        PETSC_ERR_ARG_NULL,
        "Output DofLayout pointer must not be null.");

    void *context = nullptr;

    PetscCall(
        DMShellGetContext(
            dm,
            &context));

    if (context == nullptr)
    {
        SETERRQ(
            PetscObjectComm(
                reinterpret_cast<PetscObject>(dm)),
            PETSC_ERR_ARG_WRONGSTATE,
            "DMShell has no CpGrid DofLayout context.");
    }

    *layout =
        static_cast<DofLayout *>(context);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
DofLayout::createGlobalVectorCallback_(
    DM dm,
    Vec *vector)
{
    PetscFunctionBegin;

    DofLayout *layout = nullptr;
    PetscCall(context_(dm, &layout));

    PetscCall(
        layout->createGlobalVectorImpl_(
            vector));

    /*
     * DMCreateGlobalVector() 要求自定义 DM 的创建回调把当前 DM
     * 显式关联到返回的 Vec。PETSc debug 构建会通过 VecGetDM()
     * 检查这一点。
     */
    PetscCall(
        VecSetDM(
            *vector,
            dm));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
DofLayout::createLocalVectorCallback_(
    DM dm,
    Vec *vector)
{
    PetscFunctionBegin;

    DofLayout *layout = nullptr;
    PetscCall(context_(dm, &layout));

    PetscCall(
        layout->createLocalVectorImpl_(
            vector));

    /*
     * DMCreateLocalVector() 同样要求返回的 Vec 关联当前 DMShell。
     */
    PetscCall(
        VecSetDM(
            *vector,
            dm));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
DofLayout::createMatrixCallback_(
    DM dm,
    Mat *matrix)
{
    PetscFunctionBegin;

    DofLayout *layout = nullptr;
    PetscCall(context_(dm, &layout));

    PetscCall(
        layout->createMatrixImpl_(
            matrix));

    /*
     * DMCreateMatrix() 在 PETSc debug 构建中会检查 Mat 是否关联 DM。
     */
    PetscCall(
        MatSetDM(
            *matrix,
            dm));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
DofLayout::createGlobalVectorImpl_(
    Vec *vector) const
{
    PetscFunctionBegin;

    PetscCheck(
        vector != nullptr,
        mesh_.communicator(),
        PETSC_ERR_ARG_NULL,
        "Output global Vec pointer must not be null.");

    PetscCall(
        VecCreateMPI(
            mesh_.communicator(),
            dofMap_.ownedDofCount(),
            dofMap_.globalDofCount(),
            vector));

    PetscCall(
        VecSetBlockSize(
            *vector,
            dofMap_.dofPerCell()));

    PetscCall(
        VecSetLocalToGlobalMapping(
            *vector,
            localToGlobalMapping_));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
DofLayout::createLocalVectorImpl_(
    Vec *vector) const
{
    PetscFunctionBegin;

    PetscCheck(
        vector != nullptr,
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_NULL,
        "Output local Vec pointer must not be null.");

    PetscCall(
        VecCreateSeq(
            PETSC_COMM_SELF,
            dofMap_.localSnapshotDofCount(),
            vector));

    /*
     * 本地快照可能以一个 cell block 为自然访问单位。
     * ghost 列表按完整 DOF 排列，因此 block size 只作为性能/显示提示。
     */
    PetscCall(
        VecSetBlockSize(
            *vector,
            dofMap_.dofPerCell()));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
DofLayout::createMatrixImpl_(
    Mat *matrix) const
{
    PetscFunctionBegin;

    PetscCheck(
        matrix != nullptr,
        mesh_.communicator(),
        PETSC_ERR_ARG_NULL,
        "Output Mat pointer must not be null.");

    if (!sparsity_)
    {
        sparsity_ =
            std::make_unique<Sparsity>(
                dofMap_);
    }

    PetscCall(
        sparsity_->createMatrix(
            matrix));

    PetscFunctionReturn(PETSC_SUCCESS);
}

void DofLayout::createMapping_(
    const std::vector<PetscInt> &globalIndexValues)
{
    PetscCallAbort(
        mesh_.communicator(),
        ISLocalToGlobalMappingCreate(
            mesh_.communicator(),
            1,
            static_cast<PetscInt>(
                globalIndexValues.size()),
            globalIndexValues.data(),
            PETSC_COPY_VALUES,
            &localToGlobalMapping_));
}

void DofLayout::createShell_()
{
    PetscCallAbort(
        mesh_.communicator(),
        DMShellCreate(
            mesh_.communicator(),
            &dm_));

    PetscCallAbort(
        mesh_.communicator(),
        DMShellSetContext(
            dm_,
            this));

    PetscCallAbort(
        mesh_.communicator(),
        DMShellSetCreateGlobalVector(
            dm_,
            createGlobalVectorCallback_));

    PetscCallAbort(
        mesh_.communicator(),
        DMShellSetCreateLocalVector(
            dm_,
            createLocalVectorCallback_));

    PetscCallAbort(
        mesh_.communicator(),
        DMShellSetCreateMatrix(
            dm_,
            createMatrixCallback_));
}

VecScatter DofLayout::createScatter_(
    bool localToGlobal,
    const std::vector<PetscInt> &globalIndexValues) const
{
    Vec globalTemplate = nullptr;
    Vec localTemplate = nullptr;
    IS globalIndices = nullptr;
    IS localIndices = nullptr;
    VecScatter scatter = nullptr;

    PetscCallAbort(
        mesh_.communicator(),
        createGlobalVectorImpl_(&globalTemplate));
    PetscCallAbort(
        PETSC_COMM_SELF,
        createLocalVectorImpl_(&localTemplate));

    std::vector<PetscInt> localIndexValues(
        globalIndexValues.size());
    std::iota(
        localIndexValues.begin(),
        localIndexValues.end(),
        static_cast<PetscInt>(0));

    PetscCallAbort(
        PETSC_COMM_SELF,
        ISCreateGeneral(
            PETSC_COMM_SELF,
            static_cast<PetscInt>(globalIndexValues.size()),
            globalIndexValues.data(),
            PETSC_COPY_VALUES,
            &globalIndices));
    PetscCallAbort(
        PETSC_COMM_SELF,
        ISCreateGeneral(
            PETSC_COMM_SELF,
            static_cast<PetscInt>(localIndexValues.size()),
            localIndexValues.data(),
            PETSC_COPY_VALUES,
            &localIndices));

    if (localToGlobal)
    {
        PetscCallAbort(
            mesh_.communicator(),
            VecScatterCreate(
                localTemplate,
                localIndices,
                globalTemplate,
                globalIndices,
                &scatter));
    }
    else
    {
        PetscCallAbort(
            mesh_.communicator(),
            VecScatterCreate(
                globalTemplate,
                globalIndices,
                localTemplate,
                localIndices,
                &scatter));
    }

    PetscCallAbort(PETSC_COMM_SELF, ISDestroy(&globalIndices));
    PetscCallAbort(PETSC_COMM_SELF, ISDestroy(&localIndices));
    PetscCallAbort(mesh_.communicator(), VecDestroy(&globalTemplate));
    PetscCallAbort(PETSC_COMM_SELF, VecDestroy(&localTemplate));

    return scatter;
}

void DofLayout::createGlobalToLocalScatter_(
    const std::vector<PetscInt> &globalIndexValues)
{
    globalToLocalScatter_ =
        createScatter_(false, globalIndexValues);
    PetscCallAbort(
        mesh_.communicator(),
        DMShellSetGlobalToLocalVecScatter(
            dm_,
            globalToLocalScatter_));
}

void DofLayout::ensureLocalToGlobalScatter_() const
{
    if (localToGlobalScatter_ != nullptr)
        return;

    const auto globalIndexValues =
        dofMap_.localToGlobalIndices();
    localToGlobalScatter_ =
        createScatter_(true, globalIndexValues);
    PetscCallAbort(
        mesh_.communicator(),
        DMShellSetLocalToGlobalVecScatter(
            dm_,
            localToGlobalScatter_));
}

Vec DofLayout::createGlobalVector() const
{
    Vec vector = nullptr;

    PetscCallAbort(
        mesh_.communicator(),
        DMCreateGlobalVector(
            dm_,
            &vector));

    return vector;
}

Vec DofLayout::createLocalVector() const
{
    Vec vector = nullptr;

    PetscCallAbort(
        mesh_.communicator(),
        DMCreateLocalVector(
            dm_,
            &vector));

    return vector;
}

Vec DofLayout::createLocalVector(
    Vec global) const
{
    if (global == nullptr)
    {
        throw std::invalid_argument(
            "createLocalVector(global) requires a valid global Vec.");
    }

    Vec local =
        createLocalVector();

    try
    {
        globalToLocal(
            global,
            local,
            INSERT_VALUES);
    }
    catch (...)
    {
        PetscCallAbort(
            PETSC_COMM_SELF,
            VecDestroy(&local));
        throw;
    }

    return local;
}

Vec DofLayout::borrowLocalVector() const
{
    Vec local = nullptr;
    PetscCallAbort(
        mesh_.communicator(),
        DMGetLocalVector(
            dm_,
            &local));
    return local;
}

Vec DofLayout::borrowLocalVector(
    Vec global) const
{
    if (global == nullptr)
    {
        throw std::invalid_argument(
            "borrowLocalVector(global) requires a valid global Vec.");
    }

    Vec local = borrowLocalVector();
    try
    {
        globalToLocal(
            global,
            local,
            INSERT_VALUES);
    }
    catch (...)
    {
        restoreLocalVector(local);
        throw;
    }
    return local;
}

void DofLayout::restoreLocalVector(
    Vec &local) const noexcept
{
    if (local == nullptr)
        return;

    PetscCallAbort(
        mesh_.communicator(),
        DMRestoreLocalVector(
            dm_,
            &local));
}

void DofLayout::globalToLocal(
    Vec global,
    Vec local,
    InsertMode mode) const
{
    if (global == nullptr ||
        local == nullptr)
    {
        throw std::invalid_argument(
            "globalToLocal requires valid Vec objects.");
    }

    PetscCallAbort(
        mesh_.communicator(),
        DMGlobalToLocalBegin(
            dm_,
            global,
            mode,
            local));

    PetscCallAbort(
        mesh_.communicator(),
        DMGlobalToLocalEnd(
            dm_,
            global,
            mode,
            local));
}

void DofLayout::localToGlobal(
    Vec local,
    Vec global,
    InsertMode mode) const
{
    if (global == nullptr ||
        local == nullptr)
    {
        throw std::invalid_argument(
            "localToGlobal requires valid Vec objects.");
    }

    ensureLocalToGlobalScatter_();

    PetscCallAbort(
        mesh_.communicator(),
        DMLocalToGlobalBegin(
            dm_,
            local,
            mode,
            global));

    PetscCallAbort(
        mesh_.communicator(),
        DMLocalToGlobalEnd(
            dm_,
            local,
            mode,
            global));
}

Mat DofLayout::createMatrix() const
{
    Mat matrix = nullptr;

    PetscCallAbort(
        mesh_.communicator(),
        DMCreateMatrix(
            dm_,
            &matrix));

    return matrix;
}

void DofLayout::release_() noexcept
{
    /*
     * DMShell 内部保存 scatter context 的引用，先销毁 DM，再释放本类引用。
     */
    if (dm_ != nullptr)
    {
        PetscCallAbort(
            mesh_.communicator(),
            DMDestroy(&dm_));
    }

    if (globalToLocalScatter_ != nullptr)
    {
        PetscCallAbort(
            mesh_.communicator(),
            VecScatterDestroy(
                &globalToLocalScatter_));
    }

    if (localToGlobalScatter_ != nullptr)
    {
        PetscCallAbort(
            mesh_.communicator(),
            VecScatterDestroy(
                &localToGlobalScatter_));
    }

    if (localToGlobalMapping_ != nullptr)
    {
        PetscCallAbort(
            mesh_.communicator(),
            ISLocalToGlobalMappingDestroy(
                &localToGlobalMapping_));
    }
}

} // namespace MPMC
