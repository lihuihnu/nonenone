# Natural 模块说明书

## 1. 职责

`models/include/natural` 是储层多相多组分物理核心。它把编译期 `Indices`、`FluidSystem`、相态管理、物性、蓄积、通量、井源和局部方程组合成 grid-independent cell kernel；PETSc 适配独立放在 `natural/petsc`。

## 2. 外部 API

纯 C++ 推荐入口：

```cpp
#include <natural/natural.hpp>
```

主要公共类型：

- `FluidSystem<Indices>`：流体模型与 constitutive callbacks。
- `CompositionalMixture` / `CubicEquationOfState`。
- `CubicThreePhaseFlash<Indices>`：PR/SW/CPA 共用的三相 flash；`PengRobinsonThreePhaseFlash` 仅保留为兼容别名。
- `PhasePresence`、`CompositionalPhase`。
- `NaturalCellKernel<Indices>`：cell-level 物理入口。

PETSc 推荐入口：

```cpp
#include <natural/petsc/natural_petsc.hpp>
```

## 3. 内部子模块

- `thermo`：EOS/flash。
- `state`：primary/secondary state 与 phase management。
- `properties`：密度/黏度等物性。
- `physics`：accumulation/flux/well/特殊物理。
- `assembly`：局部 residual。
- `kernel`：上述部分的统一组合器。
- `petsc`：分布式 runtime 和 Jacobian assembly。

## 4. 整体逻辑

```text
PrimaryVariables
  -> CellStateCodec
  -> phase equilibrium / phase state
  -> CellPropertyEvaluator
  -> CellProperties
  -> NaturalCellKernel
       -> accumulation
       -> face flux
       -> well source
       -> algebraic equilibrium/closure equations
  -> CellResidualBlock
```

## 5. 扩展规则

新增物理时先判断它属于 thermodynamics、constitutive property、conservation term 还是 runtime concern。不要在 `NaturalPetscRuntime` 中直接写新的 EOS 或 constitutive formula。

## 6. 不可破坏约束

- 纯 Natural 层保持 PETSc/grid independent。
- 一套物理行为只保留一个生产实现。
- Fully-compositional 三相和 Legacy 路径必须通过 `if constexpr` 保持各自语义。

## 7. 验证

`natural_core_test`、`natural_extended_test`、`natural_conservation_test`、三相/AD/井相关 unit tests。
