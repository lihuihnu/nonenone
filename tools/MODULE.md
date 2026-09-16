# Tools 模块说明书

## 1. 职责

`tools` 是 MPMC_SCW 的**离线扩展工具层**。它复用已经验证的生产物理内核，但不参与 reservoir time loop、PETSc SNES 或网格离散。适合放置相图、PVT 扫描、参数研究、结果后处理、输入转换等“需要项目物理能力、但不是求解器本体”的功能。

当前功能中心是 phase-diagram / flash-plot 工具：用与储层模拟器相同的 `CubicThreePhaseFlash` 批量扫描 P/T/z，输出 MATLAB 友好的 CSV，并生成配套 MATLAB 绘图脚本。除基础相区图外，当前版本提供连续泡点/露点包络、critical-locus estimate、unrestricted O/G/W 相出现边界、phase-quality contour、P/T flash profile，并可对 PR/SW/CPA 三个 backend 生成对比图。

依赖方向：

```text
 tools/example
      -> tools
      -> indices + models/natural/thermo
      -> common/ad (由生产热力学间接使用)

 production modules  -X->  tools
```

生产模块不得反向 include `tools`。

## 2. 外部 API

公共头文件位于 `tools/include/tools/`。

### `phase_diagram.hpp`

- `AxisSpacing::{Linear,Logarithmic}`：扫描轴间距。
- `ScanAxis{minimum,maximum,points,spacing}`：包含两端点的一维扫描轴。
- `FlashMode::{Unrestricted,Restricted}`：是否允许生产 flash 自由创建所有相。
- `PhaseDiagramFlashPolicy`：相图采样时的 phase-set 策略与失败处理。
- `EnvelopeOptions`：泡点/露点边界提取与临界候选估计配置。
- `WaterCriticalReference`：IAPWS ordinary-water critical coordinate，仅作扫描/绘图参考。
- `MultiphaseBoundaryOptions` / `MultiphaseBoundaryKind`：O/G/W phase-onset 提取配置与类别。
- `PhaseDiagramSample<Indices>`：一个 P-T-z 点及完整 `ThreePhaseFlashResult`。
- `PressureTemperatureMap<Indices>`：规则 P-T 网格。
- `PressureCompositionMap<Indices>`：固定 T 下，在两个总体组成之间线性插值并扫描压力。
- `TernaryCompositionMap<Indices>`：固定 P/T 下的三元重心网格。
- `EnvelopePoint`：一个 bubble/dew 坐标点。
- `PressureTemperatureEnvelope<Indices>`：固定总体组成下的 bubble/dew 曲线与近似 critical point。
- `PressureCompositionEnvelope<Indices>`：固定温度下沿组成路径的 bubble/dew 曲线。
- `CriticalLocusPoint<Indices>` / `CriticalLocusPath<Indices>`：沿组成路径扫描得到的临界候选曲线，并携带 envelope-width / composition-distance residual。
- `PressureTemperaturePhaseBoundaries<Indices>`：unrestricted O/G/W P-T 相出现边界。
- `PhaseDiagramSampler<Indices>`：主采样器。
- `writeCsv(...)`：对所有 map/envelope/locus 的重载 CSV writer。

主采样接口：

```cpp
PhaseDiagramSampler<Indices> sampler(flash, policy);

auto pt  = sampler.pressureTemperature(z, temperatureAxis, pressureAxis);
auto px  = sampler.pressureComposition(T, pressureAxis, zStart, zEnd, pathAxis);
auto tri = sampler.ternaryComposition(P, T, {componentA, componentB, componentC}, subdivisions);

auto ptEnvelope = sampler.pressureTemperatureEnvelope(z, temperatureAxis, pressureAxis, options);
auto pxEnvelope = sampler.pressureCompositionEnvelope(T, pressureAxis, zStart, zEnd, pathAxis, options);
auto critical   = sampler.criticalLocus(pathAxis, zStart, zEnd, temperatureAxis, pressureAxis, options);
auto onset      = sampler.pressureTemperaturePhaseBoundaries(pt, boundaryOptions);
```

`PhaseDiagramSampler` **不实现任何新 EOS/flash**；每个点都调用生产 `CubicThreePhaseFlash`。

### `eos_parameter_regression.hpp`

离线、backend-agnostic 的有界参数回归工具。production solver 不依赖它。

- `RegressionParameter`：参数名、初值、上下界、初始步长。
- `PatternSearchOptions`：最大迭代、扩张/收缩系数、停止步长和 improvement tolerance。
- `RegressionResult`：收敛状态、目标函数调用数、初始/最终目标、拟合值与最终步长。
- `boundedPatternSearch(...)`：确定性的有界 derivative-free coordinate/pattern search。
- `writeRegressionResultCsv(...)`：输出可追溯拟合结果。

Objective 完全由调用者定义，因此可以拟合：

- PR/SW/CPA constant 或 temperature-dependent cubic BIP；
- component volume translation `c_i`（前提是目标数据包含 density/volume）；
- CPA `a0,b,c1,epsilon,beta`；
- CPA cross-association `epsilon_ij,beta_ij`。

工具不会自动把 fitted values 写回 case，也不提供统计置信区间；拟合参数必须经过独立 validation set 验证后再人工进入正式流体配置。

### `matlab_phase_diagram.hpp`

- `writePressureTemperatureMatlab(...)`
- `writePressureCompositionMatlab(...)`
- `writeTernaryMatlab(...)`
- `writePressureTemperatureEnvelopeMatlab(...)`
- `writePressureCompositionEnvelopeMatlab(...)`
- `writeCriticalLocusMatlab(...)`
- `writePressureTemperatureQualityMatlab(...)`
- `writePressureFlashProfileMatlab(...)`
- `writeTemperatureFlashProfileMatlab(...)`
- `writeMultiphaseBoundariesMatlab(...)`
- `writeBackendEnvelopeComparisonMatlab(...)`
- `writeCriticalLocusComparisonMatlab(...)`

生成的 MATLAB 脚本直接 `readtable()` 读取 CSV。三元图还使用 CSV 中的相组成绘制代表性 O-G tie lines；若存在 `phase_code=7`，会绘制一个代表性的三相 tie triangle。

## 3. CSV 公共约定

所有基础相图 CSV 都包含：

- `pressure_Pa`, `pressure_bar`, `temperature_K`
- `status_code`
- `phase_code`, `phase_count`, `iterations`
- `beta_o`, `beta_g`, `beta_w`
- `So`, `Sg`, `Sw`
- `Z_o`, `Z_g`, `Z_w`
- `rho_m_o`, `rho_m_g`, `rho_m_w`
- `z_<component>`
- `xo_<component>`, `yg_<component>`, `xw_<component>`
- `Kg_<component>=yg/xo`, `Kw_<component>=xw/xo`

相态码与 `PhasePresence::bits()` 完全一致：

```text
0 = flash failed / non-converged
1 = O
2 = G
3 = O+G
4 = W
5 = O+W
6 = G+W
7 = O+G+W
```

`status_code`：`0=converged`, `1=flash returned non-converged`, `2=exception captured while continueAfterFailure=true`。失败点不会被静默删除。

各 scan 额外增加自己的坐标列：

- P-T：`temperature_index`, `pressure_index`
- P-composition：`path_index`, `pressure_index`, `path_fraction`
- ternary：`ternary_a/b/c`, `ternary_x/y`

新增 envelope / critical CSV：

- P-T envelope：`temperature_K`, `has_dew`, `dew_pressure_*`, `has_bubble`, `bubble_pressure_*`, `envelope_width_*`, `critical_point`，并给出 `critical_relative_envelope_width`, `critical_oil_gas_composition_L1`, `critical_estimate_score`
- P-composition envelope：`path_fraction`, `temperature_K`, `has_dew`, `has_bubble`, `envelope_width_*`
- critical locus：`path_fraction`, `valid`, `critical_temperature_K`, `critical_pressure_*`, `critical_envelope_width_*` 和同一组 critical-estimate residual
- phase-onset boundary：`boundary_code/name`, `lower_phase_code`, `upper_phase_code`

## 4. 内部 API

`MPMC::tools::detail` 中的以下内容只服务 writer/sampler：

- composition normalization / interpolation
- CSV column-name sanitization
- common header/row writer
- MATLAB string quoting/preamble/scatter helper
- envelope boundary helper逻辑

这些 helper 不是跨模块 API，不应被其他生产模块直接依赖。

## 5. 整体逻辑

```text
EOS / CubicThreePhaseFlash
        |
        v
PhaseDiagramSampler
        |
        +-- scan P,T at fixed z
        +-- scan P and composition path at fixed T
        +-- scan ternary z at fixed P,T
        +-- extract bubble/dew lines from monotone pressure slices
        +-- estimate critical candidates from envelope coalescence + O/G composition distance
        +-- extract oil/gas/water presence-bit boundaries
        |
        v
PhaseDiagramSample[] / EnvelopePoint[] / CriticalLocusPoint[]
        |
        +-- writeCsv()
        |
        +-- write*Matlab()
        v
CSV + .m
```

bubble/dew 提取不是重新解闪蒸，而是先用已有样点识别相边界区间，再在区间内继续调用生产 flash 做一维 pressure bisection/refinement。

critical locus 当前是**显式标记的近似数值工具**：对路径上每个总体组成提取 fixed-z P-T envelope，再用相对 bubble/dew span 与 O/G phase-composition L1 distance 构成 residual score，选择最接近 envelope coalescence 的候选点。CSV 保留这些 residual，便于论文中筛选可信点；它不是独立的严格 Michelsen criticality-equation / Hessian solver，因此图例统一写 `critical estimate`。

## 6. Restricted 与 Unrestricted 的选择

普通烃类 VLE 图推荐：

```cpp
PhaseDiagramFlashPolicy::restricted(
    PhasePresence(PhasePresence::oilBit | PhasePresence::gasBit));
```

因为这类问题只研究 O/G 平衡，不应该把 fully-compositional 架构中的第二 liquid slot 当成一个独立水相。

H2O-CO2-hydrocarbon、PR/SW/CPA 三相图通常使用：

```cpp
PhaseDiagramFlashPolicy::unrestricted();
```

并给 `ThreePhaseFlashOptions::waterComponent` 设置真实 H2O 索引。

## 7. 扩展规则

新增 tools 功能时：

1. 能复用生产 API 就不要复制公式。
2. 通用逻辑放 `tools/include/tools/`；具体研究代码放 `tools/example/<application>/`。
3. 工具不得成为 production dependency。
4. 大批量输出优先使用 CSV；绘图软件层不进入求解器。
5. 新功能增加纯 C++ unit regression；若工具仅调用纯热力学，不应强制引入 PETSc。
6. 若加入更严格的 continuation/critical solver，应与当前 sampling/extraction API 并存，而不是破坏现有 CSV/脚本接口。

未来适合扩展：strict bubble/dew continuation、严格 critical-point tracing、PVT tables、parameter sweep、带 covariance/uncertainty 的更高级 regression、flash batch benchmark、CMG/WinProp comparison exporters。

## 8. 三相含水与三 backend 使用约定

H2O–CO2–hydrocarbon 研究应把 `ThreePhaseFlashOptions::waterComponent` 设置为真实水索引，并用 `PhaseDiagramFlashPolicy::unrestricted()` 生成真实 O/G/W map、phase-onset 和 flash profiles。

论文里常见的 bubble/dew curve 属于 conventional VLE 定义。若第三个 water-rich phase 会参与，应另外建立 `restricted(O|G)` sampler 作为 **O/G VLE projection**；不要把 water-onset curve 改名为 bubble/dew line。三个 EOS backend 均复用这一规则。

`WaterCriticalReference` 的 647.096 K / 22.064 MPa 仅用于坐标参考。扫描到 850 K / 600 bar 只代表工具数值覆盖，不自动扩大 PR/SW/CPA 参数的实验有效域。

## 9. 构建与验证

```bash
cd tools
make clean
make
make run-example
```

正式单元测试由 `test/tools_phase_diagram_test` 覆盖，另见 `tools/example/MODULE.md`。
