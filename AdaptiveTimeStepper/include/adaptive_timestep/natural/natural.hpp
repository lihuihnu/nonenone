#pragma once

/**
 * @file natural.hpp
 * @brief Natural 三相多组分模型的稳定公共入口。
 *
 * 每次尝试前保存已接受的解、相态和井状态；拒绝时恢复，仅在时间步接受后提交历史量。
 */
#include <adaptive_timestep/natural/natural_petsc_backend.hpp>
