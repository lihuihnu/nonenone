# Indices 模块说明书

## 1. 职责

`indices` 是 Natural 模型的编译期布局中心。它把“模型有哪些物理功能”转换为固定的相索引、主变量索引、方程索引、相态数组索引以及 AD 导数维数。

## 2. 外部 API

### `CompositionalModelConfig`

```cpp
template <int NumComponents, bool HasWater, bool HasWellUnknown,
          bool HasAqueousCO2Dissolution, bool HasAdsorption,
          bool HasLandTrapping, PhaseBehaviorModel PhaseBehavior>
struct CompositionalModelConfig;
```

`PhaseBehaviorModel` 当前区分：

- `LegacyOilGasWithIndependentWater`
- `FullyCompositionalThreePhase`

### `Indices<Config, UseAD>`

对外公开：

- 编译期能力：`numComponents`、`numPhases`、`hasWater`、`hasWellUnknown`、`hasAdsorption`、`hasLandTrapping`、`fullyCompositionalThreePhase` 等。
- `Indices::Phase`：oil/vapor/water 相索引。
- `Indices::Primary`：pressure、各相独立组成、饱和度、可选 BHP 的位置。
- `Indices::Equation`：component balance、fugacity equilibrium、closure、well control 的位置。
- `Indices::PhaseState`：持久 secondary-state 布局。
- `ValueType`：double 或固定维度 AD 标量。
- `ADIndices<Config>`、`ScalarIndices<Config>`：常用别名。

## 3. 内部 API

`detail::makeIndexArray()` 以及 `Primary`/`Equation` 内部的 `legacy...`、`full...` 中间常量仅用于构造布局。其他模块应使用最终公共名称，不应依赖中间 offset。

## 4. 整体逻辑

```text
CompositionalModelConfig
   -> 编译期 feature flags
   -> Indices
        -> numPrimaryVariables
        -> numEquations
        -> numPhaseStateVariables
        -> AD derivative dimension
   -> Natural / PETSc / Output 共用同一布局
```

Fully-compositional 三相模式保持 `numPrimaryVariables == numEquations`；Legacy 模式保持旧版数值布局兼容。

## 5. 扩展规则

新增物理未知量时必须同时回答：

1. 哪种 Config 才启用？
2. 主变量放在哪里？
3. 对应方程是什么？
4. phase-state 是否需要持久字段？
5. 旧配置的数值索引是否必须保持？
6. `numPrimaryVariables == numEquations` 是否仍成立？

## 6. 不可破坏约束

Legacy 布局是兼容基线；不能为了新功能随意重新编号。Fully-compositional 模式不能同时启用 legacy aqueous-CO2 独立方程。

## 7. 验证

`test/src/unit/indices_test.cpp`，并由所有 Natural/AD 编译测试间接覆盖。
