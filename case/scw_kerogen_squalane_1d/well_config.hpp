/**
 * @file well_config.hpp
 * @brief 单重质拟组分超临界水驱算例的井配置入口。
 */
#pragma once

#include "../scw_kerogen_common/benchmark_common.hpp"

namespace WellConfig
{
inline static constexpr const auto &wells = ScwKerogen1D::wells;
} // namespace WellConfig
