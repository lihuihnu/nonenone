#pragma once

/**
 * @file adaptive_timestep.hpp
 * @brief 自适应隐式时间步模块的公共入口，汇总策略、后端和井控制接口。
 *
 * 模块本身只负责步长策略和统计；具体后端负责非线性求解、状态回滚/提交以及可选的井控制更新。
 */
#include <adaptive_timestep/core/config.hpp>
#include <adaptive_timestep/core/events.hpp>
#include <adaptive_timestep/core/exceptions.hpp>
#include <adaptive_timestep/core/nonlinear_stagnation.hpp>
#include <adaptive_timestep/core/policy.hpp>
#include <adaptive_timestep/core/solve_result.hpp>
#include <adaptive_timestep/core/statistics.hpp>
#include <adaptive_timestep/core/stepper.hpp>
