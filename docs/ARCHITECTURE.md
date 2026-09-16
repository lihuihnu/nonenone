# MPMC_SCW architecture

## 1. Design rule

项目保持单向依赖和单一生产实现：物理模型不依赖网格，网格不实现 EOS，PETSc 只承担并行向量、SNES/KSP 和装配驱动；离线 `tools` 复用 production 热力学，不反向进入 simulator。

```text
indices
  ↓
ad + common
  ↓
models/natural ───── models/well
  ↓                   ↓
models/natural/petsc ← grid
  ↓
AdaptiveTimeStepper
  ↓
output

production thermo ──→ tools (offline only)
```

## 2. Modules

### `indices/`

定义模型能力、主变量、方程和相状态索引。公共入口：

```cpp
#include <indices/indices.hpp>
```

核心类型：`CompositionalModelConfig`, `ADIndices`。

### `ad/`

固定维度自动微分基础类型。不能依赖 PETSc 或具体物理模型。

### `common/`

数学、插值、单位、控制台输出、PETSc 通用 I/O 等共享能力。禁止放入 EOS、井或具体网格业务逻辑。

### `models/natural/`

纯 Natural 物理核心，主要子模块：

```text
thermo/       EOS, stability, flash
properties/   density, viscosity, relperm, capillary, LBC...
physics/      accumulation/flux/source physical terms
state/        secondary state and accepted-step state
kernel/       grid-independent cell/face physics entry
assembly/     local residual/Jacobian assembly helpers
petsc/        PETSc runtime adapter
```

公共入口：

```cpp
#include <natural/natural.hpp>
#include <natural/petsc/natural_petsc.hpp>
```

`NaturalCellKernel` 是纯单元物理入口；`NaturalCpGridRuntime` / `NaturalStructuredGridRuntime` 才接 PETSc 和网格。

当前相态提交采用两级事务语义：`state/` 先在单 cell 内保证 failed Flash 不改变传入 primary/phase-state；PETSc runtime 再把 owned-cell 候选保存在 scratch 中，通过 MPI 全局判定后一次性提交。可恢复热力学失败不会生成部分更新的分布式状态，而是在下一次 SNES residual 评价中作为 function-domain failure 交给 nonlinear/timestep 层。Newton trial state 的 EOS/物性评价异常采用同一 SNES domain 语义：当前 residual 清零并结束该次 nonlinear solve，外层时间步事务负责回滚和缩步。

StructuredGrid 与 CpGrid 的高频短期 local Vec 均使用各自 DM 的 borrow/restore 临时向量池；长期状态 Vec 仍保持显式 create/destroy ownership。

v67.10 起，`natural_petsc_runtime.hpp` 只保留稳定模板类壳和数据成员，具体实现按 diagnostics/state/residual/cache/wells/jacobian/scaling/history 放在 `natural/petsc/detail/*.inc`。这些 `.inc` 不是公共入口，只用于降低单文件职责和维护风险。three-phase Flash 同样保留 `CubicThreePhaseFlash` 公共类壳，core/nonlinear/restricted/allocation/canonicalization/stability 实现位于 `natural/thermo/detail/*.inc`。

### `models/well/`

井定义、Peaceman WI、Rate/BHP 控制、schedule 和控制切换。公共入口：

```cpp
#include <well/well.hpp>
```

### `grid/`

- `StructuredGrid`：规则结构网格；
- `CpGrid`：corner-point 网格、ghost/owned DOF、岩石属性；
- GRDECL：支持 `SPECGRID/DIMENS`, `COORD`, `ZCORN`, `ACTNUM`, `PORO`, `PERMX/Y/Z`, `INCLUDE` 的当前子集。

网格层不允许包含 EOS、相平衡或井控逻辑。

### `AdaptiveTimeStepper/`

只负责时间推进策略、尝试/接受/回滚和井控重解协调。Natural 后端负责提交/恢复模型状态。

### `output/`

统一 CSV、井历史、质量库存、逐组分质量守恒和可选 VTK。结果文件只在运行时产生。

### `case/`

正式超算算例。每个算例至少：

```text
<case>.cpp
case_config.hpp
well_config.hpp
```

可选 `case_fluid.hpp` 和 `README.md`。`case/Makefile` 自动发现，不维护中央算例列表。

井的物理参数仍由各 case 的 `well_config.hpp` 定义；`case/include/case/well_factory.hpp` 只负责把这些参数统一翻译成 `WellSpecification`、处理 StructuredGrid completion/perforation 样板和可选井 limits/schedule。Panfili 等特殊时间调度仍留在对应 case。

### `test/`

唯一测试入口。业务模块不再放独立 test Makefile。

### `tools/`

PETSc-free 离线扩展：相图、P/T/z 扫描、EOS 参数回归和 MATLAB 脚本生成。每个 thermodynamic sample 都调用 production `CubicThreePhaseFlash`。

## 3. Recommended public headers

| 需求 | 入口 |
|---|---|
| 模型配置/索引 | `indices/indices.hpp` |
| Natural 核心 | `natural/natural.hpp` |
| Natural + PETSc | `natural/petsc/natural_petsc.hpp` |
| 井 | `well/well.hpp` |
| 时间步 | `adaptive_timestep/adaptive_timestep.hpp` |
| Natural 时间步后端 | `adaptive_timestep/natural/natural.hpp` |
| 输出 | `output/output.hpp` |
| PETSc/CpGrid 输出 | `output/petsc/output_petsc.hpp` |
| 相图工具 | `tools/phase_diagram.hpp` |

不要从算例直接依赖 `natural/assembly/*`、`natural/petsc/*_backend.hpp` 等内部文件，除非正在开发对应模块。

## 4. Fully compositional variable layout

在 `FullyCompositionalThreePhase` 下，H2O 与其它组分统一参与 EOS。典型 unknown layout 包含：

- pressure；
- 两个独立 saturation，第三个由闭合得到；
- O/G/W 各相的 `N-1` 独立摩尔分数；
- 可选井未知量。

约束包括逐组分质量守恒、相间化学势/逸度平衡和体积/组成闭合。具体索引由 `ADIndices<Config>` 生成，调用者不应硬编码位置。

## 5. Main runtime path

```text
case_config / well_config
        ↓
Grid + rock loading
        ↓
FluidSystem / thermodynamic backend
        ↓
P-T-z initialization / accepted state
        ↓
Natural runtime installs SNES residual/Jacobian
        ↓
AdaptiveTimeStepper
   ├─ nonlinear solve
   ├─ well-control iteration
   ├─ accept / rollback
   └─ output / diagnostics
```

## 6. Extension rules

1. 一个物理行为只保留一个 production 实现。
2. 禁用物理使用 `if constexpr`，避免空 runtime adapter。
3. 固定组分数路径优先 `std::array`，减少动态分配。
4. PETSc `VecGetArray*` 通过 RAII view 成对管理。
5. 新 EOS 通过 thermodynamic backend 接口接入，不改 grid、well 或 output。
6. 新工具只允许依赖 production physics；production 不得 include `tools`。
7. 新算例只放 `case/`；新回归只放 `test/`。
8. 新公共 API 必须同步目标模块 `MODULE.md`。
