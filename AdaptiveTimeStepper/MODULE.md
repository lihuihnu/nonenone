# AdaptiveTimeStepper 模块说明书

## 1. 职责

该模块实现与具体模型无关的事务式自适应隐式时间推进：尝试、非线性求解、井控重解、接受/拒绝、rollback、dt growth/cut 和统计。

## 2. 外部 API

纯核心推荐：

```cpp
#include <adaptive_timestep/adaptive_timestep.hpp>
```

主要类型：

- `AdaptiveTimeStepConfig`：固定输出步、最小 dt、cut/growth、easy/difficult iteration thresholds、最大 retries、最大井控循环。
- `AdaptiveTimeStepPolicy`：给定 remaining time 产生 attempt plan，并根据结果更新 dt。
- `AdaptiveTimeStepper<Backend,Observer>`：主控制器。
- `NonlinearSolveResult`、`NewtonIterationRecord`、`AdaptiveAttemptPlan`。
- `AdaptiveTimeStepStatistics`。
- 事件：`AdaptiveAttemptEvent`、`AdaptiveNonlinearSolveEvent`、`AdaptiveAcceptedEvent`、`AdaptiveRejectedEvent`。

Natural/PETSc bridge：

- `NaturalAdaptiveBackend<Runtime,...>`。
- `PetscSnesDriver`。
- `WellControlCycle` / `NoWellControlCycle`。

## 3. 内部 API

`AdaptiveTimeStepper::advanceTo_()`、retry logic、policy 的 dt clipping/growth 判断属于内部实现。Backend 的契约比具体类型重要：必须支持 begin/solve/update controls/accept/reject/current time 等事务操作。

## 4. 整体逻辑

```text
fixed output target
  -> policy.plan(remaining)
  -> backend.beginAttempt()
  -> backend.solve()
  -> well control settled?
       no -> switch + re-solve same dt
  -> converged?
       no -> reject + rollback + cut dt + retry
       yes -> accept + commit history + choose next dt
  -> 到达 fixed output time
```

## 5. 扩展规则

新增模型只需实现 backend contract；不要复制 stepper。新增 timestep policy 应保持 accepted/rejected 事务语义和 fixed-output alignment。

## 6. 不可破坏约束

- rejected attempt 不得提交 phase/well/Land/history。
- 不能跳过固定 output target。
- 井控切换必须在同一物理 timestep 内重解直至 settled 或达到上限。

## 7. 验证

`adaptive_timestep_test.cpp`、`adaptive_well_test.cpp`。
