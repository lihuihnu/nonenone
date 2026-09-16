/**
 * @file case_fluid.hpp
 * @brief Sun-2024 Exp.12 派生算例相对渗透率入口。
 */
#pragma once

#include "benchmark_common.hpp"

namespace BenchmarkCaseFluid
{
template <class Indices>
void apply(MPMC::FluidSystem<Indices> &fluid)
{
    Sun2024Exp12::applyRelativePermeability(fluid);
}
} // namespace BenchmarkCaseFluid
