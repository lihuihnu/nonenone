#pragma once

#include <common/console.hpp>

/**
 * @file output.hpp
 * @brief 网格无关输出模块公共入口。
 *
 * PETSc/CpGrid 数据采集适配器只从 `output/petsc` 头文件暴露，核心格式化/写文件
 * 代码保持 PETSc 无关，便于单元测试和其他后端复用。
 */
#include <output/core/format.hpp>
#include <output/core/writer.hpp>
#include <output/core/text_file.hpp>
#include <output/core/types.hpp>
#include <output/metrics/component_totals.hpp>
#include <output/metrics/component_mass_balance.hpp>
#include <output/well/well_history.hpp>

