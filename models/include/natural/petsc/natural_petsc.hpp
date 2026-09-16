/**
 * @file natural_petsc.hpp
 * @brief Natural/PETSc 模块公共入口。
 */
#pragma once

// Natural + PETSc 公共入口。用户通常只需要：
// 1. 选择 CpGrid 或 StructuredGrid runtime；
// 2. 用 installNaturalCallbacks() 绑定 SNES。
#include <natural/petsc/callbacks.hpp>
#include <natural/petsc/natural_cpgrid_runtime.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
