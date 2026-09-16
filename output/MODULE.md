# Output 模块说明书

## 1. 职责

`output` 负责把求解器状态转换为可追踪结果和诊断文件。它只读 simulation state，不参与 Newton 更新或 phase switching。

## 2. 外部 API

### Core

推荐入口：`output/output.hpp`。

- `OutputWriter`：solution/phase-state/mass-series CSV 基础 writer。
- `OutputTextFile`：RAII 文本文件。
- `OutputCellIdSpace`、`OutputCellSelection`、`MassSeries`。

### Metrics

- `ComponentTotals`
- `ComponentMassBalanceSnapshot`
- `ComponentMassBalanceLedger`

### Simulation

推荐入口：`output/simulation_output.hpp`。

- `ModelInventoryOutput`
- `ComponentMassBalanceOutput`
- `ReservoirDiagnosticsOutput`
- `DetailedWellOutput`
- `WellControlSwitchOutput`

### PETSc/CpGrid

- `CpGridOutputAccess`
- `CpGridSaver`
- `OutputPetscOwnedReadView`

### VTK

- `CpGridVtkWriter`、`VtkCellOrdering`。

## 3. 内部 API

CSV header 生成、sample row formatting、MPI/global data collection、file-open state 属于内部细节。Case 应通过 `CaseOutput` 或公开 Output 类型，不应直接操作 writer 内部 stream。

## 4. 整体逻辑

```text
accepted solution / runtime diagnostics
  -> input-order collection
  -> inventory / well / mass-balance / reservoir metrics
  -> CSV rows
  -> optional VTK
```

质量守恒 ledger 的累计注入/产出只在 accepted time progression 更新，避免 rejected attempts 重复积分。

## 5. 扩展规则

新增输出应先明确是“原始状态字段”“派生诊断”还是“时间积分 ledger”。派生量尽量从 runtime 的统一 global query 获取，避免再次做昂贵、重复且可能不一致的物性计算。

## 6. 不可破坏约束

- 输出不得修改 solution/phaseState。
- CpGrid 默认输出保持 input/canonical ordering。
- CSV 列名和单位一旦作为正式结果接口使用，修改必须明确记录兼容影响。

## 7. 验证

`output_test.cpp`、`output_feature_test.cpp`、`output_well_test.cpp`，以及正式 case 的 CSV 对比。
