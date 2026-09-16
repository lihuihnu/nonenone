# Case 模块说明书

## 1. 职责

`case` 是工程装配层：选择 Config/Indices、建立 FluidSystem、建立 Grid、配置 wells、建立 Natural runtime/SNES、启动 AdaptiveTimeStepper 和 Output。它不应复制底层物理公式。

## 2. 标准算例外部接口

每个正式算例目录保持：

```text
case/<name>/
  <name>.cpp
  case_config.hpp
  well_config.hpp
```

可附加 reference CSV、case-specific fluid helper 和 README。

公共装配 helper 位于 `case/include/case/`；`well_config_schema.hpp` 统一井配置数据结构，各 case 只保留参数与兼容别名。机器环境与通用运行命令见 `docs/HPC_RUN.md` 和 `case/README.md`。

### `case_config.hpp`

作为默认参数唯一来源，定义：Grid、Fluid、Numerics、Time、Output、Land/Adsorption/Dissolution 等静态配置。

### `well_config.hpp`

建立标准 `WellSpecification` 列表或 schedule updater。

### 共用 API

`case/include/case`：

- `readRunOptions<Config>()`
- `makeRuntimeOptions<Indices,Runtime,Config>()`
- `makeTimeStepConfig<Config>()`
- `NaturalSolver<Runtime>`：SNES/Vec/Mat RAII。
- `CaseOutput<...>`：统一输出组合器。
- `runTimeLoop(...)`：标准时间推进入口。
- `printSimulationConfiguration(...)` / `printRunSummary(...)`。
- `makeFluidSystem`/thermodynamic selector helpers。
- `SolverDiagnosticsObserver`。

## 3. Thermodynamic empirical parameter 配置

`makeFluidSystem<Indices,Config>()` 保持旧 `Config::Fluid::binaryInteraction` 常数矩阵为默认。若需要更丰富经验参数，可按以下方式 **opt-in**。

### 3.1 常数 BIP（旧接口，继续支持）

```cpp
inline static constexpr std::array<std::array<double,N>,N>
binaryInteraction = {/* ... */};
```

### 3.2 线性温度相关 BIP

```cpp
inline static constexpr double binaryInteractionReferenceTemperature = 350.0; // K
inline static constexpr std::array<std::array<double,N>,N>
binaryInteractionTemperatureSlope = {/* 1/K */};
```

factory 自动构造：

\[
k_{ij}(T)=k_{ij}(T_{ref})+s_{ij}(T-T_{ref}),
\]

其中 `binaryInteraction[i][j]` 就是 `Tref` 处的值，`s_ij` 单位为 `1/K`。

### 3.3 任意温度相关 BIP callback

若 case 定义：

```cpp
static double cubicBinaryInteractionCoefficient(
    int i, int j, double temperature);
```

则 callback 优先于线性 convenience 配置和常数矩阵。可在里面实现论文/实验拟合的多项式、倒温度或分段相关式。返回值必须有限。

当前接口只依赖 `T`，没有 `P/x` 参数，也不提供非等温 AD 的 `dk/dT`。

### 3.4 Component volume translation

```cpp
inline static constexpr std::array<double,N>
componentVolumeTranslation = {/* m3/mol */};
```

本项目符号约定：

\[
v=v_{EOS}-\sum_i x_i c_i.
\]

只影响 density/phase volume/saturation/flow-volume quantities，不修改 `Z`、fugacity 或 phase equilibrium。默认不声明该字段等价于全部 `c_i=0`。

SW 还可提供仅作用于 Water 热力学角色的附加平移：

```cpp
inline static constexpr std::array<double,N>
soreideWhitsonAqueousVolumeTranslation = {/* m3/mol */};
```

工厂仅在 `thermodynamicModel=SoreideWhitson` 时读取该字段。通用平移与水相附加
平移在线性摩尔体积修正中相加；油相和气相不读取水相附加项。该接口适合表达
Garcia 溶解 CO2 表观摩尔体积等密度闭合，但仍不得用于修改 SW 逸度或 BIP。

### 3.5 优先级与兼容性

```text
arbitrary cubicBinaryInteractionCoefficient(i,j,T)
    > linear slope + Tref
    > legacy constant binaryInteraction
```

SW water-rich H2O-solute pair仍由已有 `soreideWhitsonAqueousWaterBip(component,T,salinity)` 控制；generic BIP 只作为非水相/非水 pair 的 cubic coefficient。CPA 的 cubic SRK contribution 同样消费 generic BIP，而 association/cross-association 参数接口保持原样。

## 4. 内部 API

CSV column helper、console formatting、PETSc option temporary override、case validation traits 属于装配层内部。算例主函数不应直接重复这些逻辑。

## 5. 整体逻辑

```text
main
  -> PetscInitialize
  -> read RunOptions
  -> build Grid + report
  -> build FluidSystem
  -> build Runtime options
  -> NaturalPetscRuntime
  -> wells / schedule
  -> initial solution PTz or explicit state
  -> initializeHistory
  -> NaturalSolver
  -> runTimeLoop
  -> final CSV + summary
```

`case_config.hpp` 提供默认值；PETSc options 只做运行时临时覆盖。

## 6. 新增算例

1. 新建标准三文件。
2. 在 Config 里选择 `CompositionalModelConfig` 和 EOS。
3. 使用已有 fluid/grid/runtime factory。
4. 仅在 case-specific reference 数据确实需要时增加额外文件。
5. `case/Makefile` 自动发现目录，不维护中央 case 列表。

## 7. 不可破坏约束

- 普通 `make` 不覆盖用户已有 `run.sh`。
- case 不直接实现 EOS/flash/flux/井控。
- 默认参数只在 `case_config.hpp` 保留一份。

## 8. 验证

代表性正式算例：`3p6c_original`、`DQcase`、`3p4c_pr_reservoir`、`ma2021_5c_three_phase_reservoir`、`panfili2025_case3_fullphysics`。
