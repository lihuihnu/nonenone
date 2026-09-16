/**
 * @file petsc_owned_read_view.hpp
 * @brief PETSc owned 区只读 Vec 的共享 RAII 实现。
 */
#pragma once

#include <petscvec.h>

#include <stdexcept>

namespace MPMC::petsc::detail
{

/**
 * @brief 统一管理 owned range、只读数组获取与归还，不做 ghost scatter 或重排。
 */
class PetscOwnedReadViewCore
{
protected:
    explicit PetscOwnedReadViewCore(Vec vector, const char *invalidVectorMessage)
        : vector_(vector)
    {
        if (vector_ == nullptr)
            throw std::invalid_argument(invalidVectorMessage);

        const MPI_Comm comm = PetscObjectComm(
            reinterpret_cast<PetscObject>(vector_));
        PetscCallAbort(comm, VecGetOwnershipRange(vector_, &begin_, &end_));
        PetscCallAbort(comm, VecGetArrayRead(vector_, &array_));
    }

    PetscOwnedReadViewCore(const PetscOwnedReadViewCore &) = delete;
    PetscOwnedReadViewCore &operator=(const PetscOwnedReadViewCore &) = delete;
    PetscOwnedReadViewCore(PetscOwnedReadViewCore &&) = delete;
    PetscOwnedReadViewCore &operator=(PetscOwnedReadViewCore &&) = delete;

    ~PetscOwnedReadViewCore() noexcept
    {
        if (vector_ != nullptr && array_ != nullptr)
        {
            PetscCallAbort(
                PetscObjectComm(reinterpret_cast<PetscObject>(vector_)),
                VecRestoreArrayRead(vector_, &array_));
        }
    }

public:
    [[nodiscard]] const PetscScalar *data() const noexcept { return array_; }
    [[nodiscard]] PetscInt begin() const noexcept { return begin_; }
    [[nodiscard]] PetscInt end() const noexcept { return end_; }
    [[nodiscard]] PetscInt localSize() const noexcept { return end_ - begin_; }

private:
    Vec vector_{nullptr};
    const PetscScalar *array_{nullptr};
    PetscInt begin_{0};
    PetscInt end_{0};
};

} // namespace MPMC::petsc::detail
