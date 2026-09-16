# Natural Kernel 子模块说明书

## 1. 职责

`NaturalCellKernel<Indices>` 是纯 C++ Natural 物理的 facade。它组合 state decode、property evaluation、face flux 和 residual assembly，使 PETSc runtime 不需要了解各物理文件之间的细节。

## 2. 外部 API

主要方法类别：

- 从 primary values + phase state 构造 `CellState`。
- 计算 `CellProperties`。
- 计算 face mass flux。
- 组装一个 cell 的 `CellResidualBlock`。

具体数据类型来自 `state/cell_state.hpp`、`assembly/cell_residual.hpp`。

## 3. 内部 API

`CellStateCodec`、`CellPropertyEvaluator` 等成员对象是 kernel 内部协作者。上层 runtime 不应重复执行相同物理步骤。

## 4. 整体逻辑

```text
NaturalCellKernel
  ├─ decode state
  ├─ evaluate properties
  ├─ evaluate face physics
  └─ assemble residual
```

## 5. 扩展规则

如果一个新功能是“每个 cell 必须统一执行”的物理步骤，应先在对应 physics/state 层实现，再由 kernel 组合，而不是直接把长公式写进 kernel。

## 6. 不可破坏约束

Kernel 不能持有 PETSc Vec/Mat；输入输出必须是普通 C++ 数据结构/AD 标量。

## 7. 验证

`natural_core_test.cpp`、`natural_extended_test.cpp`、`well_natural_test.cpp`。
