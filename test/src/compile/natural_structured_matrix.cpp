/**
 * @file natural_structured_matrix.cpp
 * @brief 编译矩阵测试：验证 `natural_structured_matrix` 所覆盖模板/后端组合能够实例化。
 */
#include <indices/indices.hpp>
#include <natural/fluid_system.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>

#include <array>

/**
 * @brief 编译期强制实例化 StructuredGrid backend 的功能配置矩阵。
 *
 * 与 CpGrid 配置矩阵保持一致：基础、溶解、Land、吸附、无井以及全组分油气水三相。
 */
template <class Config>
void forceNaturalStructuredGridConfiguration(
    MPMC::StructuredGridCore &grid,
    const MPMC::FluidSystem<
        MPMC::ADIndices<Config>> &fluid,
    Vec solution,
    Vec residual,
    Mat jacobian)
{
    using Indices = MPMC::ADIndices<Config>;
    using Grid = MPMC::StructuredGridCore;
    using Runtime =
        MPMC::NaturalStructuredGridRuntime<
            Indices,
            Grid>;

    typename Runtime::Options options{};

    if constexpr (Indices::hasAdsorption)
    {
        options.standardGasDensity =
            std::array<double, Indices::numComponents>{};
        options.standardGasDensity->fill(1.0);
    }

    Runtime runtime(grid, fluid, options);
    if constexpr (Indices::fullyCompositionalThreePhase)
    {
        std::array<double, Indices::numComponents> z{};
        z.fill(1.0 / static_cast<double>(Indices::numComponents));
        runtime.initializeUniformFromPTZ(solution, 1.0e7, fluid.temperature, z);
    }
    else
    {
        runtime.initializePhaseStateFromSolution(solution);
    }
    runtime.updateState(solution);
    runtime.initializeHistory(solution);
    (void)runtime.formFunction(solution, residual);
    (void)runtime.formJacobian(solution, jacobian, jacobian);
    runtime.limitNewtonStep(solution, residual);
}

using Base =
    MPMC::CompositionalModelConfig<
        6, true, true,
        false, false, false>;
using Dissolution =
    MPMC::CompositionalModelConfig<
        6, true, true,
        true, false, false>;
using Land =
    MPMC::CompositionalModelConfig<
        6, true, true,
        false, false, true>;
using Adsorption =
    MPMC::CompositionalModelConfig<
        6, true, true,
        false, true, false>;
using NoWell =
    MPMC::CompositionalModelConfig<
        6, true, false,
        false, false, false>;
using FullyCompositional =
    MPMC::CompositionalModelConfig<
        4, true, true,
        false, false, false,
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;

#define MPMC_FORCE_STRUCTURED_CONFIG(ConfigType)                           \
    template void forceNaturalStructuredGridConfiguration<ConfigType>(     \
        MPMC::StructuredGridCore &,                                        \
        const MPMC::FluidSystem<                                           \
            MPMC::ADIndices<ConfigType>> &,                                \
        Vec, Vec, Mat)

MPMC_FORCE_STRUCTURED_CONFIG(Base);
MPMC_FORCE_STRUCTURED_CONFIG(Dissolution);
MPMC_FORCE_STRUCTURED_CONFIG(Land);
MPMC_FORCE_STRUCTURED_CONFIG(Adsorption);
MPMC_FORCE_STRUCTURED_CONFIG(NoWell);
MPMC_FORCE_STRUCTURED_CONFIG(FullyCompositional);

#undef MPMC_FORCE_STRUCTURED_CONFIG
