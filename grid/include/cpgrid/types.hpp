/**
 * @file types.hpp
 * @brief CpGrid 公共类型、索引和轻量别名。
 */
#pragma once

#include <mpi.h>

namespace MPMC
{

/**
 * @brief CpGrid 使用的 MPI 通信器类型。
 */
using CpCommunicator = MPI_Comm;

} // namespace MPMC
