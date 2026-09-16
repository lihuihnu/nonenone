#pragma once

/**
 * @file well.hpp
 * @brief 井模型公共入口，汇总定义、控制、调度和 Peaceman 井指数。
 *
 * Flow sign convention is uniform across this module: injection is positive and
 * production is negative. Reservoir coupling is implemented by Natural's well source.
 */
#include <well/control.hpp>
#include <well/manager.hpp>
#include <well/peaceman.hpp>
#include <well/perforation.hpp>
#include <well/schedule.hpp>
#include <well/specification.hpp>
#include <well/state.hpp>
#include <well/types.hpp>
#include <well/vertical_well.hpp>
