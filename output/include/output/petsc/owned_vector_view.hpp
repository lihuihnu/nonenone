/**
 * @file owned_vector_view.hpp
 * @brief PETSc 全局向量 owned 单元块的只读视图。
 */
#pragma once

#include <common/petsc_owned_read_view.hpp>

#include <petscvec.h>

namespace MPMC
{

class OutputPetscOwnedReadView final
    : private petsc::detail::PetscOwnedReadViewCore
{
public:
    explicit OutputPetscOwnedReadView(Vec vector)
        : PetscOwnedReadViewCore(
              vector,
              "OutputPetscOwnedReadView requires a valid PETSc Vec.")
    {
    }

    OutputPetscOwnedReadView(const OutputPetscOwnedReadView &) = delete;
    OutputPetscOwnedReadView &operator=(const OutputPetscOwnedReadView &) = delete;
    OutputPetscOwnedReadView(OutputPetscOwnedReadView &&) = delete;
    OutputPetscOwnedReadView &operator=(OutputPetscOwnedReadView &&) = delete;
    ~OutputPetscOwnedReadView() noexcept = default;

    using PetscOwnedReadViewCore::data;
    using PetscOwnedReadViewCore::begin;
    using PetscOwnedReadViewCore::end;
    using PetscOwnedReadViewCore::localSize;
};

} // namespace MPMC
