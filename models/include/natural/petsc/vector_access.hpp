/**
 * @file vector_access.hpp
 * @brief PETSc Vec 的安全数组映射与只读/可写访问辅助。
 */
#pragma once

#include <common/petsc_owned_read_view.hpp>

#include <petscvec.h>

#include <stdexcept>

namespace MPMC
{

/**
 * @brief PETSc global Vec 当前 rank owned 区的只读 RAII 视图。
 *
 * 本类不做 ghost scatter；它只包装 `VecGetOwnershipRange + VecGetArrayRead`，
 * 用于 updateState/history 等只访问
 * 本 rank owned DOF 的代码。析构保证 Restore，即使中间 C++ 物理计算抛异常。
 */
class PetscOwnedReadView final
    : private petsc::detail::PetscOwnedReadViewCore
{
public:
    explicit PetscOwnedReadView(Vec vector)
        : PetscOwnedReadViewCore(
              vector,
              "PetscOwnedReadView requires a valid Vec.")
    {
    }

    PetscOwnedReadView(const PetscOwnedReadView &) = delete;
    PetscOwnedReadView &operator=(const PetscOwnedReadView &) = delete;
    PetscOwnedReadView(PetscOwnedReadView &&) = delete;
    PetscOwnedReadView &operator=(PetscOwnedReadView &&) = delete;
    ~PetscOwnedReadView() noexcept = default;

    using PetscOwnedReadViewCore::data;
    using PetscOwnedReadViewCore::begin;
    using PetscOwnedReadViewCore::end;
};

/**
 * @brief PETSc global Vec 当前 rank owned 区的可写 RAII 视图。
 */
class PetscOwnedWriteView final
{
public:
    explicit PetscOwnedWriteView(Vec vector)
        : vector_(vector)
    {
        if (vector_ == nullptr)
            throw std::invalid_argument(
                "PetscOwnedWriteView requires a valid Vec.");

        [[maybe_unused]] const MPI_Comm comm =
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    vector_));

        PetscCallAbort(
            comm,
            VecGetOwnershipRange(
                vector_,
                &begin_,
                &end_));

        PetscCallAbort(
            comm,
            VecGetArray(
                vector_,
                &array_));
    }

    PetscOwnedWriteView(const PetscOwnedWriteView &) = delete;
    PetscOwnedWriteView &operator=(const PetscOwnedWriteView &) = delete;
    PetscOwnedWriteView(PetscOwnedWriteView &&) = delete;
    PetscOwnedWriteView &operator=(PetscOwnedWriteView &&) = delete;

    ~PetscOwnedWriteView() noexcept
    {
        if (vector_ != nullptr && array_ != nullptr)
        {
            PetscCallAbort(
                PetscObjectComm(
                    reinterpret_cast<PetscObject>(
                        vector_)),
                VecRestoreArray(
                    vector_,
                    &array_));
        }
    }

    [[nodiscard]] PetscScalar *data() noexcept
    {
        return array_;
    }

    [[nodiscard]] const PetscScalar *data() const noexcept
    {
        return array_;
    }

    /** @brief 在对象析构前显式归还 Vec 数组，便于随后调用其它 Vec 操作。 */
    void restore()
    {
        if (vector_ != nullptr && array_ != nullptr)
        {
            PetscCallAbort(
                PetscObjectComm(reinterpret_cast<PetscObject>(vector_)),
                VecRestoreArray(vector_, &array_));
        }
    }

    [[nodiscard]] PetscInt begin() const noexcept
    {
        return begin_;
    }

    [[nodiscard]] PetscInt end() const noexcept
    {
        return end_;
    }

private:
    Vec vector_{nullptr};
    PetscScalar *array_{nullptr};
    PetscInt begin_{0};
    PetscInt end_{0};
};


/**
 * @brief 临时把 Vec 整体乘以 -1，并在析构时恢复原符号。
 *
 * 遵循 Natural `updateSol` 的符号约定；即使 limiter 抛异常也不会把
 * 调用方的 Newton step 永久留在反号状态。
 */
class PetscVectorSignFlipGuard final
{
public:
    explicit PetscVectorSignFlipGuard(Vec vector)
        : vector_(vector)
    {
        if (vector_ == nullptr)
            throw std::invalid_argument(
                "PetscVectorSignFlipGuard requires a valid Vec.");

        PetscCallAbort(
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    vector_)),
            VecScale(
                vector_,
                -1.0));
    }

    PetscVectorSignFlipGuard(const PetscVectorSignFlipGuard &) = delete;
    PetscVectorSignFlipGuard &operator=(const PetscVectorSignFlipGuard &) = delete;
    PetscVectorSignFlipGuard(PetscVectorSignFlipGuard &&) = delete;
    PetscVectorSignFlipGuard &operator=(PetscVectorSignFlipGuard &&) = delete;

    ~PetscVectorSignFlipGuard() noexcept
    {
        if (vector_ != nullptr)
        {
            PetscCallAbort(
                PetscObjectComm(
                    reinterpret_cast<PetscObject>(
                        vector_)),
                VecScale(
                    vector_,
                    -1.0));
        }
    }

private:
    Vec vector_{nullptr};
};

} // namespace MPMC
