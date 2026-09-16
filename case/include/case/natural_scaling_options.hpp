/**
 * @file natural_scaling_options.hpp
 * @brief Natural 方程缩放的 PETSc 命令行覆盖与日志输出。
 */
#pragma once

#include <petscsys.h>

#include <type_traits>

namespace MPMC::cases
{

namespace detail
{
template <class Numerics, class = void>
struct HasNaturalScalingCaseDefaults : std::false_type {};

template <class Numerics>
struct HasNaturalScalingCaseDefaults<Numerics, std::void_t<
    decltype(Numerics::pressureScale),
    decltype(Numerics::compositionScale),
    decltype(Numerics::saturationScale),
    decltype(Numerics::massResidualScale),
    decltype(Numerics::fugacityResidualScale),
    decltype(Numerics::closureResidualScale),
    decltype(Numerics::rateWellResidualFloor)>> : std::true_type {};
} // namespace detail

/** Apply optional case-level scales before PETSc command-line overrides. */
template <class Numerics, class RuntimeOptions>
void applyNaturalScalingCaseDefaults(RuntimeOptions &runtimeOptions)
{
    if constexpr (detail::HasNaturalScalingCaseDefaults<Numerics>::value)
    {
        runtimeOptions.scaling.pressureScale = Numerics::pressureScale;
        runtimeOptions.scaling.compositionScale = Numerics::compositionScale;
        runtimeOptions.scaling.saturationScale = Numerics::saturationScale;
        runtimeOptions.scaling.massResidualScale = Numerics::massResidualScale;
        runtimeOptions.scaling.fugacityResidualScale = Numerics::fugacityResidualScale;
        runtimeOptions.scaling.closureResidualScale = Numerics::closureResidualScale;
        runtimeOptions.scaling.rateWellResidualFloor = Numerics::rateWellResidualFloor;
    }
}

/**
 * @brief 按现有命令行约定覆盖 Natural runtime scaling 配置。
 *
 * 保持原算例的选项读取顺序、默认值、validate() 调用和 [SCALE] 输出格式。
 * 这里只改变数值坐标缩放参数，不改变任何物理方程或方程零点。
 */
template <class RuntimeOptions>
void applyNaturalScalingPetscOptions(RuntimeOptions &runtimeOptions, PetscMPIInt rank)
{
    PetscBool scalingEnabled = PETSC_TRUE;
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetBool(
            nullptr, nullptr, "-natural_scaling", &scalingEnabled, nullptr));
    runtimeOptions.scaling.enabled = (scalingEnabled == PETSC_TRUE);

    PetscReal pressureScale = runtimeOptions.scaling.pressureScale;
    PetscReal compositionScale = runtimeOptions.scaling.compositionScale;
    PetscReal saturationScale = runtimeOptions.scaling.saturationScale;
    PetscReal massScale = runtimeOptions.scaling.massResidualScale;
    PetscReal fugacityScale = runtimeOptions.scaling.fugacityResidualScale;
    PetscReal closureScale = runtimeOptions.scaling.closureResidualScale;
    PetscReal rateFloor = runtimeOptions.scaling.rateWellResidualFloor;

    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetReal(
            nullptr, nullptr, "-natural_pressure_scale", &pressureScale, nullptr));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetReal(
            nullptr, nullptr, "-natural_composition_scale", &compositionScale, nullptr));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetReal(
            nullptr, nullptr, "-natural_saturation_scale", &saturationScale, nullptr));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetReal(
            nullptr, nullptr, "-natural_mass_residual_scale", &massScale, nullptr));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetReal(
            nullptr, nullptr, "-natural_fugacity_residual_scale", &fugacityScale, nullptr));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetReal(
            nullptr, nullptr, "-natural_closure_residual_scale", &closureScale, nullptr));
    PetscCallAbort(
        PETSC_COMM_WORLD,
        PetscOptionsGetReal(
            nullptr, nullptr, "-natural_rate_residual_floor", &rateFloor, nullptr));

    runtimeOptions.scaling.pressureScale = pressureScale;
    runtimeOptions.scaling.compositionScale = compositionScale;
    runtimeOptions.scaling.saturationScale = saturationScale;
    runtimeOptions.scaling.massResidualScale = massScale;
    runtimeOptions.scaling.fugacityResidualScale = fugacityScale;
    runtimeOptions.scaling.closureResidualScale = closureScale;
    runtimeOptions.scaling.rateWellResidualFloor = rateFloor;
    runtimeOptions.scaling.validate();

    if (rank == 0)
        PetscPrintf(
            PETSC_COMM_SELF,
            "[SCALE] enabled=%s U(p/x/S)=%.6e/%.6e/%.6e "
            "R(mass/fug/closure)=%.6e/%.6e/%.6e rate_floor=%.6e\n",
            runtimeOptions.scaling.enabled ? "ON" : "OFF",
            runtimeOptions.scaling.pressureScale,
            runtimeOptions.scaling.compositionScale,
            runtimeOptions.scaling.saturationScale,
            runtimeOptions.scaling.massResidualScale,
            runtimeOptions.scaling.fugacityResidualScale,
            runtimeOptions.scaling.closureResidualScale,
            runtimeOptions.scaling.rateWellResidualFloor);
}

} // namespace MPMC::cases
