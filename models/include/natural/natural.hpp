#pragma once

/**
 * @file natural.hpp
 * @brief Natural 三相多组分模型的稳定公共入口。
 *
 * 不依赖 PETSc 的调用方只需包含本头文件即可构造流体系统和单元计算内核；
 * PETSc 运行时接口集中在 `natural/petsc/natural_petsc.hpp`。
 */
// Natural 纯 C++ 公共入口：模型配置在 indices 模块中定义，
// 本头文件只暴露建立流体模型和单元离散所需的稳定类型。
#include <natural/compositional_mixture.hpp>
#include <natural/fluid_system.hpp>
#include <natural/kernel/cell_kernel.hpp>
#include <natural/phase_state.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <natural/thermo/phase_behavior.hpp>
#include <natural/thermo/three_phase_flash.hpp>
