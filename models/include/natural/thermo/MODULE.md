# Natural Thermodynamics 子模块说明书

## 1. 职责

`thermo` 实现等温多组分相平衡：普通 PR、Søreide–Whitson、CPA、相稳定性、二/三相 Rachford–Rice 与 PTz flash。三种 cubic backend 共享统一的经验参数扩展接口：温度相关 cubic BIP 与可选 volume translation；这些扩展默认关闭。

## 2. 外部 API

### 模型选择

- `CubicThermodynamicModel::{PengRobinson,SoreideWhitson,CubicPlusAssociation}`。
- `CompositionalPhase::{Oil,Vapor,Water}`。
- `PhasePresence`：1–7 phase mask。

### `CubicEquationOfState<Indices>`

主要对外能力：

- `mixingParameters(...)`
- `phaseMixing(...)`
- `liquidRoot(...)` / `vaporRoot(...)`
- `phaseResult(...)`
- `compressibility(...)`
- `fugacity...`
- CPA/SW option 配置
- `configureBinaryInteractionFunction(...)`
- `binaryInteractionCoefficient(i,j,T)`
- `configureVolumeTranslation(...)`
- `molarDensity(P,T,x,Z)`
- `translatedMolarVolume(P,T,x,Z)`
- `cubicPlusAssociationOptions()`：离线诊断/回归的只读 CPA 参数视图。

### Cubic BIP 的三层输入

所有 PR/SW/CPA cubic contribution 都最终调用同一个 effective coefficient：

```text
1. arbitrary callback k_ij(T)             highest precedence
2. case factory linear correlation        optional convenience layer
3. legacy binaryInteraction constant      fallback / old behavior
```

EOS 层直接支持：

```cpp
eos.configureBinaryInteractionFunction(
    [](int i, int j, double temperature) {
        return /* fitted k_ij(T) */;
    });
```

若未配置 callback，`binaryInteractionCoefficient(i,j,T)` 严格回退到 `CompositionalMixture::binaryInteractionCoefficient(i,j)`。

当前 callback 的温度是 `double`，因此它适合本项目当前等温 reservoir model 和离线 PVT/phase-diagram 扫描；它不提供能量方程所需的 `dk/dT` AD 导数。

### Søreide–Whitson BIP precedence

SW water-rich phase 保留已有 aqueous H2O-solute correlation/callback。其优先级为：

```text
water-rich H2O-i pair
  -> SW aqueousWaterBip_i(T, salinity)

non-water pair / non-aqueous phase
  -> generic cubic k_ij(T)
  -> legacy constant BIP fallback
```

因此新增的 generic `k_ij(T)` 不会覆盖已有 SW 水相专用经验相关式。

CO2-H2O 提供两个可显式选择的公开关联式：

- `co2AqueousBip(...)`：Søreide–Whitson 1992 原始参数化，保留兼容性；
- `co2AqueousBipChabab2019(...)`：Chabab et al. 2019 改进参数化，用于需要改善高 NaCl molality CO2 溶解度的算例。

两者不会隐式切换；算例必须在 `aqueousWaterBip` callback 中明确选择。

### Volume translation

可选组分常数：

\[
c_i\quad [\mathrm{m^3/mol}].
\]

采用本项目明确的符号约定：

\[
v_{\mathrm{translated}}
=\frac{ZRT}{P}-\sum_i x_i c_i,
\qquad
\rho_m=\frac{1}{v_{\mathrm{translated}}}.
\]

Volume translation 是 post-EOS density/volume correction：

- **不修改** `Z`；
- **不修改** fugacity coefficient；
- **不修改** fugacity equality / stability / flash phase split；
- 修改 molar density、mass density、由 phase moles 转 saturation 的体积换算，以及依赖 density 的 accumulation/flux/viscosity correction。

默认 `c_i=0` 时走显式 legacy fast path，避免为了数学等价重排浮点运算而改变旧参考输出。

### CPA 可拟合参数

CPA 继续支持已有完整参数：

- pure-component `a0`, `b`, `c1`；
- association `epsilon`, `beta`；
- donor/acceptor site counts；
- cubic `k_ij(T)`；
- explicit cross-association `epsilon_ij`, `beta_ij`；
- simplified / Carnahan–Starling RDF。

默认 CPA 参数与 CR-1 fallback 保持独立；cubic BIP 可随温度变化，参数拟合只通过 `tools` 离线完成。

### `CubicThreePhaseFlash<Indices>`

这是当前 canonical 三相 flash 类型，通过 `CubicEquationOfState` 统一支持 PR、SW、CPA。`PengRobinsonThreePhaseFlash<Indices>` 仅作为兼容 alias。

主要输入/输出：

- `ThreePhaseFlashOptions`
- `ThreePhaseFlashResult`
- `ThreePhaseStabilityResult`
- full PTz flash、restricted flash、stability evaluation、saturation-from-moles 更新。

当 volume translation 开启时，flash 的 phase equilibrium 仍使用未平移 EOS fugacity；仅最终 moles-to-saturation 换算使用 translated molar density。SW 可另配只作用于 Water thermodynamic role 的 `aqueousVolumeTranslation`；该附加项同样只进入密度/相体积路径，并由三相 flash、secondary-state 与流动物性使用一致的相角色。

公开 flash 返回前统一 canonicalize O/G/W public slots。普通 PR/CPA 的 O/G 方向使用 Wilson volatility；SW 内部则把 `Aqueous / HydrocarbonLiquid / Vapor` thermodynamic roles 与 public slots 解耦，单次 nonlinear solve 中 role 固定，收敛后才重排输出。若 Water 角色启用了专用物性闭包，`configureAqueousCompositionDomain()` 还要求候选组成位于该闭包的显式适用域；“水是最大摩尔组分”本身不足以把高含水烃富相移入 Water 槽。canonicalization 不允许通过改变标签偷偷改变已求解方程。

## 3. 内部 API

`three_phase_flash.hpp` 中的 `Rr2/Rr3`、log-K state、stability seeds/residual、cubic parameter cache、canonicalization helper 均为数值实现细节。外部模块不应依赖它们。

CPA 内部 association site fractions、radial distribution、density root helper 也属于实现细节。

`binaryInteractionFunction_`、`volumeTranslation_` 及 SW aqueous translation 是 EOS 实例状态；调用方应通过公开配置/查询 API 使用，不直接依赖内部存储。

## 4. 整体逻辑

```text
P,T,z
  -> normalize / initial K seeds
  -> stability trials
  -> choose phase set
  -> 2-phase or 3-phase RR/log-K nonlinear solve
  -> EOS fugacity updates
       -> effective cubic k_ij(T)
       -> SW aqueous BIP where applicable
       -> CPA association where applicable
  -> canonical oil/water role assignment
  -> final oil/gas role assignment by Wilson volatility
  -> phase compositions + beta + Z
  -> EOS molarDensity()
       -> legacy P/(RTZ), or
       -> translated 1/(ZRT/P - sum x_i c_i)
  -> saturation conversion
```

CPA 路径在 cubic contribution 外增加 Wertheim association contribution；SW 路径使用 phase-aware aqueous correlations/BIP。

## 5. 扩展规则

新增 thermodynamic backend 应尽可能继续实现 `CubicEquationOfState` 可消费的统一接口，使 flash/state/flow 无需复制。若数值算法只适用于某 EOS，应通过显式 capability/strategy 表达，而不是在上层通过 case 名称分支。

经验参数扩展遵守：

1. 常数旧参数必须继续作为 fallback；
2. 新相关式通过显式 opt-in 配置；
3. 参数拟合属于 `tools`，production EOS 只负责消费参数；
4. 如果未来增加非等温能量方程，必须重新设计 BIP callback 以提供 T 的 AD/导数语义，而不是继续把 `double T` 当作充分接口。

## 6. 不可破坏约束

- EOS 参数、BIP、association 参数、volume-shift 的单位/来源必须可追溯。
- 不通过放松 stability/flash tolerance 换性能或拟合效果。
- root selection、phase canonicalization 和 trace-component 保护是求解语义的一部分。
- volume translation 不得悄悄进入 fugacity/equilibrium equations；若未来引入不同 translated-EOS 理论，必须作为独立策略明确命名。
- 未显式启用经验参数扩展时，case 必须保持常数参数默认路径。

## 7. 验证

核心回归：`three_phase_flash_test`、`phase_role_canonicalization_test`、`sw_flash_recovery_test`、`thermodynamic_failure_modes_test` 和 `thermo_empirical_parameters_test`。

参数工具回归：`eos_parameter_regression_test`。项目级热力学设计与经验参数边界见 `../../../../docs/THERMODYNAMICS.md`。

### 当前 Flash 数值语义

PR/SW/CPA 都先执行标准 active-set/direct flash；失败后可进入共用 stable reduced-set enumeration，枚举 O+G/O+W/G+W，且只接受 restricted flash 收敛、缺失相稳定性 `valid && stable` 的候选，再按 Gibbs 选优。SW 在此之后才允许进入 exact-mass allocation/Gibbs/chemical-potential recovery、固定 thermodynamic-role permutation 和 pressure continuation。所有成功出口统一检查相分率、总体组分物料闭合与 resolved fugacity closure；stability 求解无效不能作为稳定证书。
