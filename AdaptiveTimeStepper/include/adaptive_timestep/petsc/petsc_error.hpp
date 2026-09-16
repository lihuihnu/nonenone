/**
 * @file petsc_error.hpp
 * @brief PETSc 错误码到 C++ 异常的统一转换辅助。
 */
#pragma once

#include <petscsys.h>

#include <sstream>
#include <stdexcept>
#include <string>

namespace MPMC
{

inline void throwOnPetscError(PetscErrorCode code, const char *operation)
{
    if (code == PETSC_SUCCESS)
        return;

    std::ostringstream message;
    message << operation << " failed with PETSc error code " << code << '.';
    throw std::runtime_error(message.str());
}

} // namespace MPMC
