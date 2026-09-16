# Thermodynamics and three-phase flash

## 1. Scope

当前 fully-compositional 路径使用统一的 O/G/W 三相接口，可选择：

- ordinary Peng–Robinson (PR)；
- Søreide–Whitson (SW)；
- Cubic-Plus-Association (CPA)：默认 SRK physical term，算例可显式 opt-in PR physical term。

上层 flow、well、state 和 PETSc 接口不因 backend 改变。

## 2. P-T-z initialization

给定压力、温度和总体摩尔分数 `z`：

1. 归一化总体组成；
2. 执行相稳定性测试；
3. 根据 active phase set 求单相、两相或三相平衡；
4. 得到相组成 `x^O`, `y^G`, `x^W` 与相摩尔分率 `β_O, β_G, β_W`；
5. 由 EOS 摩尔密度把相摩尔分率转换为 saturation；
6. 写入 runtime 初态。

三相基本约束：

```text
β_O + β_G + β_W = 1
z_i = β_O x_i^O + β_G y_i^G + β_W x_i^W
f_i^O = f_i^G = f_i^W       (resolved components)
```

痕量组分低于数值 composition floor 时不作为 resolved fugacity-closure 分量，但仍保持整体组成/质量约束。

## 3. Peng–Robinson

普通 PR 使用同一套 component parameters 与非水相 BIP 规则。混合吸引项：

```text
a_ij = (1 - k_ij) sqrt(a_i a_j)
```

相标签本身不切换另一套 EOS 参数，因此公开 O/G/W slot 主要是输出/物理角色表示。

## 4. Søreide–Whitson

SW 保留 PR cubic 主体，但针对 water/brine 使用专用 water alpha 和 aqueous interaction rules。工程中明确区分：

```text
non-aqueous BIP matrix
    → hydrocarbon liquid / vapor

aqueous SW BIP matrix
    → aqueous thermodynamic role
```

CO2-H2O 水相 BIP 保留 1992 原式，并新增 Chabab et al. 2019（DOI
`10.1016/j.ijggc.2019.102825`）公开改进式。二者由算例 callback 显式选择，
不会根据盐度自动切换；这既保持历史结果可复现，也避免高盐实验修正静默改变其他气体或非水相参数。

### 4.1 Thermodynamic role 与 public slot 解耦

从 v41 起，内部 thermodynamic role 为：

```text
HydrocarbonLiquid
Vapor
AqueousLiquid
```

公开存储仍为：

```text
Oil / Gas / Water
```

单次非线性求解期间 role 固定；aqueous BIP 与 liquid/vapor root 由固定 role 决定，不能因为 O/G/W 输出槽位重排而改变。收敛并验证后才一次性 canonicalize 到标准 public slots。

这消除了旧实现中的闭环：

```text
slot label → aqueous/non-aqueous parameters → fugacity → composition → slot label
```

### 4.2 Current flash fallback pipeline

当前三 EOS 都遵守“正常路径优先、失败后才恢复”的 fail-only 语义：

```text
normal active-set flash
        ↓ fail
historical direct three-phase fallback
        ↓ fail
stable reduced-set enumeration (O+G / O+W / G+W)
        ↓ only SW continues
exact-mass allocation + Gibbs basin search
+ chemical-potential LM/Newton
        ↓ fail
fixed thermodynamic-role permutation search
        ↓ fail
pressure continuation
```

其中 stable reduced-set enumeration 是 PR/SW/CPA 共用的热力学 active-set fallback：候选必须 restricted flash 收敛、缺失相稳定性 `valid && stable`，再按 dimensionless Gibbs 选优。SW 才继续进入 allocation、thermodynamic-role permutation 与 pressure continuation 等专用恢复。任何 `stability.valid=false` 的约化状态都不能作为已获得稳定性证书的成功结果返回。

v56 对 **simulator 明确指定的 SW restricted 两相 active-set** 增加边界终检：若 safeguarded two-phase 流程留下 `converged=false`，但返回状态同时满足相分率闭合、总体组成闭合和逐组分逸度平衡，则接受该物理解。若单个非水相的 liquid/vapor cubic 根已经合并，restricted G+W/O+W 保留调用者明确指定的固定 thermodynamic role；普通 unrestricted SW flash 和内部 fail-only recovery 仍使用历史 public O/G canonicalization。

### 4.3 Pressure continuation

v42 的最终 pressure continuation 只解决 globalization basin gap，不改变热力学方程。算法：

1. 目标点完整 SW 路径失败；
2. 在邻近压力搜索可收敛锚点；
3. 把公开 O/G/W 解逆映射到固定内部 role；
4. 由锚点组成构造 `K_g`, `K_w`；
5. 沿 `log(P)` 小步推进；
6. 每一步仍调用同一个三相平衡方程；
7. 到达目标后执行同样的物料、角色和 fugacity 验证。

最后的 300 K / 206.597966 bar 点用 200 bar 已收敛三相状态作 seed 时，原三相方程可以直接收敛，证明该问题属于 attraction basin，而不是 EOS 无解。

## 5. CPA

CPA backend 使用可选的 SRK 或 PR cubic reference 加 Wertheim association Helmholtz contribution；默认保持 SRK，H2O–CO2–nC10 二维 benchmark 显式选择 PR-CPA。当前主要能力：

- 4C water association scheme；
- association site fraction equations；
- cross-association / CR-1 类型规则；
- sCPA/Carnahan–Starling radial distribution；
- association pressure 与 fugacity contribution；
- 非多项式 liquid/vapor density roots；
- AD 通过隐式导数传递。

CPA association 由组成、密度和 association 参数连续控制，不通过 public `Water` slot 离散切换参数矩阵。

v57 对等温 CPA 增加两类不改变物理方程的性能路径：只依赖温度的纯组分/交叉参数可缓存；当恰有一个对称 4C 缔合组分时，site fraction 的质量作用方程直接使用解析正根，通用 cross-association 情形仍保留原 fixed-point。可选 profiler 统计 mixing/cache、association、density-root 和完整 phase evaluation 的累计调用次数与耗时。

## 6. Empirical parameter extensions

三个 backend 可使用统一的 temperature-dependent cubic BIP callback。还支持 component volume translation；CPA 参数回归由 `tools/eos_parameter_regression.hpp` 离线完成。

这些 fitted parameters 不应由工具自动写回正式 case。任何新参数都需要独立 validation set 后再人工进入 `case_config.hpp`。

## 7. Phase disappearance and reappearance

flow solve 中，当某 active phase 的 saturation 进入不可接受区域时，状态可退化到对应 reduced phase set。存在 missing phase 时重新执行 stability test；若缺失相不稳定，则进入包含新相的 flash。

有效 phase presence 共 7 种：

```text
O, G, W, O+G, O+W, G+W, O+G+W
```

状态转换由相稳定性和物理解决定，不以 public slot 名称直接决定 EOS 参数。

正常相边界仍使用窄 active-set hysteresis；若整个 SNES 已失败，v57 允许一次 fail-only rescue 在更宽探测带内尝试 restricted 相集并对同一 dt 重求。该路径仍受 missing-phase stability 保护，只有边界级数值中性状态才可保留约化相集，明显失稳必须恢复缺失相。

## 8. Density and flow coupling

EOS 给出每一 active phase 的压缩因子/摩尔密度。相摩尔分率通过相摩尔体积换算为 saturation，随后进入 mobility、Darcy flux、accumulation 和 well source。

启用锁定的 IAPWS+Garcia 富水体积闭合时，纯水基准按 IF97 区域选择：压缩液态水
使用 Region 1；`623.15<T<=863.15 K` 且不高于 B23 边界的低密度超临界水使用
Region 2；高于 B23 边界的致密超临界水使用 Region 3 Helmholtz 方程并反算密度。
当前 Sun-2024 Exp.12 派生算例的 `673.15 K, 24 MPa` 位于 Region 2。该算例
把超临界水作为独立守恒/流动相，其压力相关密度直接调用同一 IF97 区域分派器；
若 Newton 状态进入 B23 边界以上，则切换 Region 3，而不是外推 Region 2。

因此 thermodynamics 与 saturation 的关系是：

```text
P,T,z
 → flash
 → phase compositions + beta
 → EOS molar density
 → phase volume
 → saturation
 → relperm/mobility/flow
```

## 9. Non-negotiable invariants

- public slot reorder 不能改变单次 nonlinear solve 中的 SW thermodynamic role；
- fail-only recovery 不能改写原本已收敛状态；
- 新 globalization 方法不能修改 EOS/BIP 或通过放宽 residual threshold 伪造收敛；
- 所有返回的 active phase 必须满足组分物料闭合和 resolved fugacity closure；
- PR/CPA 的既有成功路径不能因 SW 修复而变化。
