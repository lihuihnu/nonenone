# Natural State 子模块说明书

## 1. 职责

`state` 管理“求解主变量”和“热力学/物性 secondary state”之间的转换，并实现相稳定、相出现/消失、restricted flash 和 Newton 更新限制。

## 2. 外部 API

- `CellState<Indices,Scalar>`：当前单元物理状态。
- `CellProperties<Indices,Scalar>`：当前单元物性缓存数据。
- `PhaseStateData<Indices>`：跨 Newton/时间步持久化的 K、beta、Z、phase presence 等 secondary state。
- `CellStateCodec<Indices>`：primary/phase-state -> `CellState`。
- `CellPropertyEvaluator<Indices>`：`CellState` -> `CellProperties`。
- `PhaseEquilibriumManager<Indices>`：Legacy oil/gas 稳定性与 flash 管理。
- `FullyCompositionalThreePhaseEquilibrium<Indices,FlashBackend>`：O/G/W 相态管理。
- `newton_limiter.hpp`：Newton update 限制。

## 3. 内部 API

`phase_detail` / `three_phase_detail` 的 normalize、trial、mask、candidate selection 等 helper 属于内部算法。不要从 PETSc/runtime 直接调用这些 helper。

## 4. Fully-compositional 主逻辑

```text
primary + previous phaseState
  -> 检查负饱和度
  -> 删除候选消失相
  -> restricted flash on active phases
  -> missing-phase stability test
       stable   -> 保持当前 active set
       unstable -> full PTz reflash
  -> 写回 compositions / K / beta / Z / phase mask
```

相消失不会改变未知量维数；assembly 使用 inactive equation replacement。

## 5. 扩展规则

新增相态策略必须把“状态判定”和“方程定义”分开：本层只决定 phase presence/secondary state，assembly 决定对应 residual。任何 phase switching 都必须可从相同输入确定性重现。

## 6. 不可破坏约束

- 不允许在检测负饱和度前 clamp 掉相消失信号。
- rejected timestep 不提交历史 phase state。
- trace component 不能通过无界 log-fugacity residual 破坏 Newton。

## 7. 验证

`three_phase_flash_test`、`three_phase_ad_test`、reservoir preflight、thermodynamic audit regressions。
