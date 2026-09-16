# Natural/PETSc 子模块说明书

## 1. 职责

该层是“纯 Natural 物理”和“分布式 PETSc 网格/求解器”的唯一主要桥梁。它负责 Vec/Mat、owned/ghost 访问、cell cache、well runtime、residual/Jacobian assembly、历史提交和全局诊断，但不重新实现 EOS/flux/well 公式。

## 2. 外部 API

推荐入口：

```cpp
#include <natural/petsc/natural_petsc.hpp>
```

主要类型：

- `NaturalPetscRuntime<Indices,Backend>`：核心 runtime。
- `NaturalStructuredGridRuntime<Indices,...>` / `NaturalCpGridRuntime<Indices,...>`：常用 alias。
- `NaturalRuntimeOptions`：dt、fugacity scaling、Land/adsorption/dissolution runtime 参数。
- `NaturalStructuredGridBackend` / `NaturalCpGridBackend`：网格差异适配。
- `NaturalGlobalWellComponentRates`、`NaturalGlobalInventory`、`NaturalGlobalDiagnostics`：MPI 全局统计。
- `callbacks.hpp`：SNES function/Jacobian/updateState/updateSol glue。
- `jacobian_factory.hpp`：按 backend 建立 Jacobian。

Runtime 的主要生命周期 API：

```text
configure wells/options
initialize solution / phase state / history
formFunction / formJacobian
updateState / updateSol
begin/reject/commit timestep
well/global inventory/diagnostics queries
input-ordered output copy
```

## 3. 内部 API

- `vector_access.hpp` 的 RAII local/owned Vec view。
- `phase_state_codec.hpp` 的 PETSc phase-state 编解码。
- `property_conversion.hpp` 的 scalar/AD conversion。
- `well_runtime.hpp` 的 runtime well helper。
- `CellCacheEntry`、static assembly cache、PETSc object-state cache guards、well source cache。

这些内部缓存不得被 case/output 直接修改。

## 4. Residual/Jacobian 主调用链

```text
PETSc solution Vec
  -> phase-state update (必要时 stability/flash)
  -> build/reuse local CellProperties cache
  -> previous accepted property cache
  -> well states/sources
  -> owned cells
       -> NaturalCellKernel residual
       -> neighbor reciprocal derivative blocks
  -> Vec assembly

same AD residual
  -> derivative blocks
  -> MatSetValuesBlocked
  -> PETSc Jacobian
```

运行时缓存静态 topology/TPFA/DOF 和 previous-state properties；current-state cache 通过 PETSc object state 判断失效，数学公式不变。

SNES residual 评价中的 Newton trial state 可能暂时离开 EOS 定义域。`naturalFormFunction` 在此时
清零残差并设置 PETSc function-domain flag，使当前 nonlinear solve 以 divergence reason 返回；
AdaptiveTimeStepper 随后回滚事务、缩小 `dt` 并重试。该恢复只作用于 residual callback 内的
trial-state 评价；初始化、配置和 Jacobian 的硬错误仍按 PETSc/C++ 错误处理。

## 5. 扩展规则

新增 runtime 功能首先判断能否放回纯 C++ 层。只有“分布式数据访问、PETSc 生命周期、MPI reduction、assembly policy”属于本层。新增 cache 必须定义明确 invalidation condition。

## 6. 不可破坏约束

- `VecGetArray*` / restore、local vector acquire/restore 必须成对。
- rejected attempt 必须恢复 solution/phaseState/well/history。
- residual trial-state domain failure 必须返回 SNES divergence，不能绕过时间步层抛出进程级错误。
- cache 不得跨状态误复用；PETSc object state 或明确 generation counter 必须参与有效性判定。
- external cell ordering 不得暴露 MPI current ordering。

## 7. 验证

纯物理可由 unit tests 覆盖；该层最终必须在目标 PETSc/MPI/METIS 环境执行 `make full`、`make run-full` 和正式 case。
