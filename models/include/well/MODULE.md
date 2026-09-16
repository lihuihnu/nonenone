# Well 模块说明书

## 1. 职责

`well` 定义与网格/求解器无关的井数据模型：井类型、控制方式、perforation、Peaceman WI、schedule、控制 limits/switching 和 converged well state。

## 2. 外部 API

推荐入口：

```cpp
#include <well/well.hpp>
```

主要类型：

- `WellType::{Injector,Producer}`。
- `WellControl`：BHP、总/相 rate 等控制类型。
- `WellPerforation`：cell、WI 等 perforation 数据。
- `WellSpecification`：井名/id、类型、控制、target、perforations、注入组成等。
- `WellState`：BHP、surface/reservoir phase rates、mass rates 等 converged state。
- `WellControlLimits`、`WellControlSelection`、`WellControlUpdate`、`WellControlSwitchReason`。
- `WellSchedule`：时间相关控制/状态更新。
- `WellManager`：按 id 管理井集合。
- `VerticalPeacemanCell` / vertical-well helper：计算 WI。

## 3. 内部 API

`control.hpp` 中 tolerance/limit comparison helper、vertical WI 计算中间量等属于内部实现。Natural coupling 由 `natural/physics/well_source.hpp` 与 `natural/petsc/well_runtime.hpp` 完成，不应反向塞进 well core。

## 4. 整体逻辑

```text
WellSpecification + perforations
  -> WI / control target
  -> Natural local properties + BHP
  -> perforation phase/component source
  -> summed WellState
  -> control limits
  -> WellControlUpdate
  -> 必要时同一 timestep 重新求解
```

统一符号：储层 source 中注入为正、生产为负；输出可同时提供非负工程 production magnitude。

## 5. 扩展规则

新增控制类型需要同步：`WellControl`、target interpretation、runtime residual、`WellState::controlledRateMagnitude`、switch logic、输出列和 tests。

## 6. 不可破坏约束

- 一个 BHP representative cell 不能被两口井共享。
- schedule/limit switch 只在明确的控制循环中修改井控制。
- well core 不直接操作 PETSc Vec。

## 7. 验证

`well_core_test.cpp`、`well_control_test.cpp`、`well_natural_test.cpp`、`adaptive_well_test.cpp`、`output_well_test.cpp`。
