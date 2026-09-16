#pragma once

/**
 * @file config.hpp
 * @brief 模块配置参数、默认值及输入合法性约束。
 */

#ifndef MESH_DIM
#define MESH_DIM 3
#endif

static_assert(
    MESH_DIM == 2 || MESH_DIM == 3,
    "CpGrid only supports MESH_DIM == 2 or MESH_DIM == 3.");
