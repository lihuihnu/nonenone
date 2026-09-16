/**
 * @file case_fluid.hpp
 * @brief H2O-CO2-nC10 二维算例的共享流体物性配置入口。
 */
#pragma once

#include "benchmark_common.hpp"

namespace BenchmarkCaseFluid
{
template <class Indices>
void apply(MPMC::FluidSystem<Indices> &fluid)
{
    BenchmarkCommon::applyCommonRelativePermeability(fluid);
}
} // namespace BenchmarkCaseFluid
