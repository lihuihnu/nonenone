# SCW–干酪根裂解产物拟组分流动算例集合

本目录用于建立“超临界水与预生成干酪根裂解产物的相平衡—黏度—多相多组分运移”数值实验集合。**当前阶段仍是零流动热力学标定；储层流动尚未获准进入。**

## 研究边界

- 主目标工况固定为 **380 °C、25 MPa**；后续保留 **360 °C、25 MPa** 作为同压力亚临界水对照候选。
- 研究对象是**预生成的干酪根裂解产物**与注入水之间的相平衡和流动耦合。
- 当前不把干酪根裂解反应、aquathermolysis、焦炭生成或能量方程纳入本算例集合。
- 拟组分边界、PR/CPA 参数、二元作用参数和黏度模型必须有可追溯实验或文献依据，不以方便计算的代表分子反向定义真实流体。

## 实验驱动的 lumping 决策

本算例不再把 `nC4 / nC10 / squalane` 或预先规定的 `C6-C14 / C15-C20 / C21+` 当作默认物理 lumping。

Zhao et al. (*Industrial & Engineering Chemistry Research*, 2023, DOI `10.1021/acs.iecr.3c02759`) 的 380 °C generated-oil Figure 6 给出直接印刷的 simulated-distillation 分布：

| lump | boiling range | recovered-oil wt% |
|---|---|---:|
| `OIL_GASOLINE` | IBP–180 °C | 0.81 |
| `OIL_DIESEL` | 180–350 °C | 23.73 |
| `OIL_MIDDLE` | 350–500 °C | 34.11 |
| `OIL_HEAVY` | >500 °C | 41.35 |

因此当前烃相 lumping topology 直接采用这四个**实验馏程区间**。详细决策见 `03_EXPERIMENT_DRIVEN_LUMPING.md`。

这组分数来自酸洗 Type-II 干酪根实验，属于 `SECONDARY_PAIRED`：它们是当前最直接的真实实验先验，但不能伪装成完整 raw-shale 主样品的同样分数。主样品如果后续获得同物理样品 simulated-distillation，应保留相同实验驱动原则并用主样品数据重新赋值。

## 当前流体拓扑

第一阶段非反应热力学模型按以下结构组织：

- `H2O`
- `OIL_GASOLINE`：IBP–180 °C
- `OIL_DIESEL`：180–350 °C
- `OIL_MIDDLE`：350–500 °C
- `OIL_HEAVY`：>500 °C

气体产物 `H2 / CO2 / CH4 / C2+` 暂不强行并入上述油相 lumps。只有在 gas + recovered-liquid 的统一质量/摩尔基准建立、且不重复计算低沸点损失后，才增加独立气体组分或 gas lump。

## 初始 pseudo-component characterization 已建立

380 °C Figure 5 的实际碳数分布被用于恢复每个实验馏程 lump 的内部 SCN 权重；Figure 6 继续作为各 lump 总质量分数的唯一直接来源。随后采用 petroleum-fraction characterization 得到第一版 EOS screening properties：

| lump | mass frac. | mole frac. | MW g/mol | Tb °C | SG | Tc K | Pc MPa | omega | Vc cm3/mol |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `OIL_GASOLINE` | 0.008100 | 0.036330 | 82.392 | 68.991 | 0.691443 | 515.231 | 3.192179 | 0.269328 | 358.490 |
| `OIL_DIESEL` | 0.237300 | 0.406770 | 215.583 | 278.493 | 0.836212 | 736.423 | 1.751155 | 0.583473 | 824.073 |
| `OIL_MIDDLE` | 0.341100 | 0.325421 | 387.349 | 430.304 | 0.896458 | 870.437 | 1.138075 | 0.896554 | 1318.853 |
| `OIL_HEAVY` | 0.413500 | 0.231479 | 660.132 | 563.696 | 0.942605 | 982.864 | 0.808714 | 1.215676 | 1793.890 |

完整推导、方法和不确定性见 `04_PSEUDOCOMPONENT_CHARACTERIZATION.md`；机器可读表为 `fluid_characterization/pseudo_component_characterization_380c.csv`，Figure 5 的 SCN characterization basis 为 `fluid_characterization/raw/figure5_380c_scn_characterization_basis.csv`。

这些性质是**实验约束 + petroleum characterization correlation 的 provisional 值**，不是直接测得的临界性质。它们可以作为零流动 PR 标定的初始纯组分参数，但不能在未经二元实验验证的情况下直接进入储层模型。

## Heavy 与 squalane 的边界

`OIL_HEAVY` 绝不默认等于 squalane。当前 >500 °C Heavy 的中心估计约为 `MW=660 g/mol`、`Tb=564 °C`、`SG=0.943`，并且实验 SARA 显示显著 resin/asphaltene character。Figure 5 可见 C38–C74 尾部也可能低估最不挥发/最极性的残余物，因此 Heavy 行明确标记为 `PROVISIONAL_HEAVY_TAIL_LOWER_BOUND_LIKE`。

squalane 仅保留为**重饱和烃 benchmark**，其参考数据独立保存在 `fluid_characterization/squalane_benchmark.csv`，不得把其纯组分性质或 H2O–squalane BIP 复制到生产 `OIL_HEAVY` 行。

## 零流动 PR 标定是进入储层前的硬门槛

在任何 60×20×1 实验尺度 slab 或其它流动算例开始前，必须分别完成：

- `H2O–OIL_GASOLINE`
- `H2O–OIL_DIESEL`
- `H2O–OIL_MIDDLE`
- `H2O–OIL_HEAVY`

的高温高压 PR 二元标定与独立验证。

每个二元体系的 `kij(T)` 只有在同时验证以下物理量后才能被接受：

1. 相数 / phase topology；
2. 两个共存相的组成；
3. 相密度或与相支相关的实验 molar-volume/PVT 信息；
4. 相界/critical-locus 位置；
5. 独立 hold-out 温压点。

只让 flash 收敛、只拟合一个富烃相端点、或者用 BIP 补偿密度错误都不能算通过。

特别是 `H2O–OIL_HEAVY`：当前 squalane 只能作为非极性重饱和烃 benchmark；atmospheric-residue 数据虽然化学上更接近 Heavy，但现有数据的压力覆盖不足以证明 25 MPa 下的 Heavy LLE/density 行为。因此 Heavy 当前明确 `BLOCKED`，**禁止为了得到两相或收敛而调 BIP**。

完整协议见 `05_ZERO_FLOW_PR_CALIBRATION.md`。实验代理数据源、published PR prior 以及储层门禁分别保存在：

- `binary_pr_calibration/source_manifest.csv`
- `binary_pr_calibration/published_pr_bip_priors.csv`
- `binary_pr_calibration/h2o_squalane_reference.csv`
- `binary_pr_calibration/reservoir_entry_gate.csv`

当前总状态：`RESERVOIR_GATE_BLOCKED`。

## PR baseline v1 已建立

当前零流动回归使用的 PR 起始参数已经固定为一套有 provenance 的 baseline，而不再使用任意 BIP seed。纯 pseudo-component 参数见 `06_PR_PARAMETER_BASELINE.md` 和 `pr_parameters/pr_pure_parameters_380c.csv`。

当前 H2O–lump optimizer seeds 为：

- Gasoline: `kref=0.5000, b=0`，仅为 Søreide–Whitson C5+ 非水相 screening prior；
- Diesel: `kref=0.6662345, b=-1274.90 K`，来自 dodecane published-PR prior；
- Middle: `kref=0.2398346, b=-477.41 K`，来自 squalane 等效碳数饱和烃 benchmark prior；
- Heavy: 保留同一 squalane 数值只作 benchmark diagnostic，`fit_enabled=0`，禁止作为生产 Heavy BIP。

所有 hydrocarbon–hydrocarbon BIP 暂取 0，直到同油样 multicomponent PVT/phase-behavior 数据证明需要非零项。

此外，Middle/Heavy 的高 `omega` 使 PR alpha-model 本身成为显著不确定性来源。PR76 仍作为统一回归基线，但最终接受前必须进行 PR78 高-omega sensitivity，避免让 `kij` 补偿纯组分 alpha 模型误差。

## CPA 独立参数体系已建立

CPA 不复用 PR 的 `a/b/alpha/kij`。当前单独维护：

- `cpa_parameters/cpa_pure_parameters_380c.csv`：CPA/SRK physical `a0/b/c1`；
- `cpa_parameters/cpa_association_scheme.csv`：4C water 与各 oil lump 的 association policy；
- `cpa_parameters/cpa_binary_matrix_screening.csv`：CPA-specific BIP；
- `cpa_parameters/cpa_entry_gate.csv`：CPA 自己的物理验收门禁；
- `cpa_parameters/source_manifest.csv`：参数来源；
- `07_CPA_PARAMETER_BASELINE.md`：完整方法与边界。

水使用仓库公开 VLLE benchmark 已验证的 Folas 4C 参数。四个油 pseudo-component 的 baseline physical term 统一使用 SRK critical mapping；Gasoline/Diesel/Middle 暂按 non-associating petroleum fractions，Heavy 也不凭 SARA 猜 association site/epsilon/beta。Heavy association 保持 `UNRESOLVED`。

CPA 的水–烃 BIP 与 PR 完全独立：Gasoline 以 water+n-hexane CPA `kij=0.044` 为 proxy prior；Diesel 暂以 water+n-decane `kij=-0.054` 为直接表格 anchor；Middle 采用 neutral prior 等待可接受数据；Heavy 为 diagnostic zero 且禁止为了 flash 收敛回归。

## 密度与黏度独立验证集已建立

相平衡标定、密度验证和黏度验证现在是三个互相隔离的证据层。验证集位于 `property_validation/`，详细规则见 `08_DENSITY_VISCOSITY_VALIDATION.md`。

- density validation：15 条已录入行，包括 3 条 IAPWS-IF97 Region-3 官方验证状态和 12 条历史未参与拟合的 squalane hold-out density 行；
- viscosity validation：34 条已录入行，包括 11 条 IAPWS-2008 官方验证点、6 条 n-hexane reference-correlation 验证点、5 条 n-undecane reference-correlation/diagnostic 点和 12 条 squalane hold-out viscosity 行；
- 所有这些行均为 `VALIDATION_ONLY` 或 implementation/diagnostic evidence，不允许用于回归 PR/CPA BIP、volume translation 或 transport 参数；
- 黏度首先使用 reference density 验证 intrinsic transport closure；只有该层通过后才允许做 EOS-density + viscosity 的 coupled validation；
- squalane 对 Middle/Heavy 仍然只是 saturated-heavy benchmark，不能让 Heavy 因此通过验证。

当前总状态仍为 `DENSITY_GATE_BLOCKED` 与 `VISCOSITY_GATE_BLOCKED`，主要缺口是 target-window 的真实油馏分/长链代理高温高压 density/viscosity 数据，尤其是 Heavy。

## 完整 0D PVT 验收已设为流动前最后一道热力学门禁

PR 与 CPA 现在必须在同一组零维状态上分别完成 PVT preflight，不能用流动结果反推哪个 EOS“合理”。验收协议见 `09_0D_PVT_ACCEPTANCE.md`，机器可读组成扫描见 `pvt_acceptance/composition_scan.csv`。

目标温度为 `360 / 374 / 380 °C`，压力为 `25–30 MPa`。组成扫描以实验四个 oil lump 的摩尔比为中心，沿 H2O overall mole fraction 从 `0.01` 扫到 `0.995`，并加入 light-enriched / heavy-enriched 两个明确标为 deterministic sensitivity 的油组成族。现在还增加了 `H2O_HEAVY_BINARY` 下维筛查族：Gasoline/Diesel/Middle 总库存严格为 0，Heavy=`1-z_H2O`，用于隔离单一 Heavy pseudo-component 的相平衡与物性行为。详细结果和物理门禁见 `16_H2O_HEAVY_BINARY_SCREENING.md`。

每个状态输出并检查：

- production flash convergence；
- final active-set stability；
- phase count / phase code；
- 两/三相组成；
- mass / molar density；
- LBC viscosity 与适用时的 IAPWS water viscosity reference；
- Oil / Gas / Water canonical phase role；
- material closure 与 active-phase fugacity closure；对总体库存为零的组分不强制无物理意义的跨相逸度等式，详见 `16_H2O_HEAVY_BINARY_SCREENING.md`。

主初始锚点固定为 `BASE oil ratio + z_H2O=0.20 + 25 MPa`，分别在 360/374/380 °C 下评估。PR 与 CPA 不要求预测相同相数；相数差异是模型结果。真正硬要求是**两个 EOS 各自都必须物理自洽**。

此外还输出 unrestricted O/G/W P–T map、oil/gas/water phase-onset envelope、restricted O/G bubble/dew projection，以及三温度下的 dense pressure-composition maps。

运行入口：

`make -C tools run-scw-kerogen-0d-pvt-acceptance`

硬门禁入口：

`make -C tools require-scw-kerogen-0d-pvt-acceptance`

即使 0D structural gate 通过，最终 flow comparison 仍必须同时满足 PR binary、CPA、density 和 viscosity 现有 gates；依赖关系见 `pvt_acceptance/flow_entry_gate.csv`。

当前严格 0D runtime gate 已经 **PASS**：

- PR registered scan：792/792；
- CPA registered scan：792/792；
- 其中新增 H2O–Heavy binary registered states：PR 198/198、CPA 198/198；
- PR/CPA target-window P–T map health：四个 family 全部 PASS；
- PR dense composition paths：10206/10206（BASE 5103 + H2O–Heavy 5103）；
- CPA dense composition paths：10206/10206（BASE 5103 + H2O–Heavy 5103）；
- 360/374/380 °C、25 MPa 的 BASE cross-EOS initial anchors：全部 PASS。

原先 CPA 在 380 °C、`z_H2O≈0.7610625`、26.00/26.25 MPa 的两个 failure states 已在**冻结 CPA 参数**的条件下定位为 active-set/continuation 数值问题并修复。详细证据见 `10_CPA_PHASE_ONSET_AUDIT.md`。这只清除了 0D 数值阻断，不代表 CPA 参数已经实验标定。

## 多孔介质已改为实验室试件定义

未来流动比较不再使用抽象“储层”岩石。新的 porous-medium contract 见 `11_LAB_POROUS_MEDIA_BASELINE.md` 与 `porous_media/`。

第一阶段优先参考 **fired Berea sandstone**：文献试件为直径 38 mm、长度 200 mm、孔隙度约 21%、水测绝对渗透率 212 mD、有效孔体积 47.6 mL；这些数值只用于定义实验尺度和采购/设计范围，真正进入模型时必须由实际试件重新测量。

二维实验仍可采用 sandstone slab 或 water-wet quartz sand-pack slab，但 slab 的 `phi / k / PV` 必须逐次装填/逐块测量，不能复制文献砂包的 2.9 D 或 13 D。

高温 SCW 条件下目前没有可接受的本体系相渗标定，因此第一阶段只允许使用 `porous_media/corey_sensitivity.csv` 中的简化 Corey **敏感性族**；任何一条都不得标成 experimental calibration。毛管压力同理，低温 Berea 的 Brooks-Corey 参数只作为 sensitivity proxy。

此外，当前 Natural face-flux 尚无显式 phase capillary-pressure closure。实验尺度流动前必须满足以下二选一：

1. 实现并验证 `Pc(S)`；
2. 用所选试件、流率、黏度和界面参数证明 `Pc=0` 是可接受的实验近似。

旧 `scw_kerogen_common` 中 `phi=0.35`、`kx=ky=1500 mD`、`kz=150 mD`、`L=1.2 m` 已明确降级为 legacy numerical values，不允许直接作为新的实验可复现实验参数。

## 井控改为实验流量尺度反推

实验 slab 的井控不再采用旧 benchmark 的“等储层体积注采”。

新的控制契约见 `13_LAB_FLOW_CONTROL_AND_PVI.md` 与：

- `porous_media/flow_control_contract.csv`
- `porous_media/darcy_flow_scale.csv`
- `porous_media/flow_control_gate.csv`

控制模式固定为：

- injector：`ReservoirTotalRate` 主控制 + `maximumBhp` 安全/有效域限制；
- producer：fixed `Bhp`。

设计压力中心取约 `28 MPa`，但注采压差不预设。实验装置/试件给出 `L / A / k / PV_eff` 且黏度验证通过后，先由

[
Q_{m inj}=r_{m PVI}PV_{m eff}
]

和

[
Delta p_{m Darcy}simrac{mu LQ}{kA}
]

反推实验可测的流速、压降和 residence time，再确定 producer BHP、expected injector pressure 与 injector maximum BHP。

当前 0D PVT 已验证约 `25–30 MPa`，因此第一阶段井压也必须保持在该窗口内；若需要超出，则必须先扩展并重跑 0D PVT gate。

最终流动结果的主时间坐标统一为**实际累计注入 PVI**：

[
mathrm{PVI}(t)=rac{int_0^t Q_{m inj,actual,res}(	au),d	au}{PV_{m eff}}.
]

如果 injector 因 maximum BHP 切换控制，PVI 使用实际注入量而不是 nominal target rate。物理时间仍保存，但不再只用“天”比较不同实验。

## 正式 SCW 温度控制组已注册

正式实验不允许只报告孤立的 380 °C 结果。主对照固定为：

- `360 °C / 28 MPa`：subcritical compressed-water control；
- `380 °C / 28 MPa`：supercritical-water test。

`374 °C / 28 MPa` 只作为 near-critical diagnostic，用于判断响应是否在临界附近出现明显非线性，不替代 360/380 主对照。

完整设计见 `14_SCW_TEMPERATURE_CONTROL_EXPERIMENT.md`。Formal pair 在 mechanism-isolation 层面锁定相同的几何、试件/均质岩石场、初始总体组成、注入组成、target PVI rate、井控、重力、Corey 与 capillary sensitivity case，只允许温度改变。

流动比较按同一 PVI 评价：

- phase topology / phase composition；
- density / viscosity / mobility；
- pressure drop / injector BHP；
- saturation/front evolution；
- component recovery、effluent enrichment 和 light/heavy selectivity。

PR 与 CPA 各自先做 `380-360` 的 within-EOS 温度差，再比较两个 EOS 对该温度效应的预测差异。

0D hard gate 现已增加 `SCW_360_380_CONTROL_PAIR_28MPA`，并输出 `scw_temperature_control_0d.csv`。最终 flow readiness 还必须通过 `porous_media/scw_temperature_control_gate.csv`。

该 28 MPa paired anchor 已实际运行并 **PASS**。当前 BASE `z_H2O=0.20` 在 PR/CPA 下 360/374/380 °C 都保持单一 Oil-role 相；360→380 °C 的初始态变化主要表现为密度和黏度降低，而不是相数突变。因此后续“SCW 效应”必须重点观察注水后水富集局部组成的 phase/mobility/selection trajectory，而不是从初始点预设一定会出现新相。

## 生产井组分升级为一级输出

正式流动不得只报告总采收率。公共 runner 现在支持独立的 \`producer_composition.csv\`，详细定义见 \`15_PRODUCER_COMPOSITION_AND_SELECTIVITY.md\`。

每口 producer 在每个固定输出状态必须给出：

- 每个 conserved component 的瞬时生产质量流率；
- 瞬时生产质量分数 \(Y_i\)；
- accepted-step 积分的累计产量 \(M_{i,\mathrm{prod}}\)；
- 以初始全局在位质量为分母的 \(RF_i\)；
- 实际累计注入 reservoir volume 与 PVI；
- 瞬时与累计 mass-basis \(E_{L/H}\)。

当前实验驱动四油 lump 的选择性诊断注册为：

\[
L=\mathrm{OIL\_GASOLINE}+\mathrm{OIL\_DIESEL},
\qquad
H=\mathrm{OIL\_HEAVY}.
\]

Middle 始终单独输出，不并入 L/H。

\[
E_{L/H}
=
\frac{Y_L/Y_H}{w_{L,0}/w_{H,0}}.
\]

这里分母使用初始**质量比**，避免把 producer 质量分数与初始摩尔分数 \(z\) 混成无物理意义的指标。

当前正式 5 组分 PVT 拓扑没有独立 gas-product lump，因此一级输出按真实 conserved components 使用 \`H2O / OIL_GASOLINE / OIL_DIESEL / OIL_MIDDLE / OIL_HEAVY\`；禁止把 Gas phase 或 Gasoline 偷换成独立 \`Gas\` component。若后续加入真实 Gas lump，必须先完成统一 gas+liquid 质量/摩尔基准。

## 当前数据审计状态

主样品仍以 Zhao et al. (*Sustainable Energy & Fuels*, 2023, DOI `10.1039/D2SE01361D`) 的完整 Chang 7 raw-shale `380 °C / 25 MPa / 4 h` 数据作为第一优先级，用于主样品油产率、SARA 和产气约束。

ACS 2023 酸洗纯干酪根数据用于提供直接的 boiling-range topology、paired mass-fraction prior 和内部 carbon-number characterization。两套样品身份保持隔离，不能把 secondary paired 数据重新标成 primary raw-shale measurement。

## 实验尺度 pseudo-3D slab（冻结）

后续流动模型保留规则 `60×20×1` 拓扑，但**不再预设物理尺寸**。完整定义见 `12_LAB_PSEUDO3D_SLAB_DESIGN.md`。

- `nx=60, ny=20, nz=1`：固定离散拓扑；
- `Lx/Ly/Lz`：由实际高温高压 slab/core-holder 的有效流动尺寸决定；
- `dx=Lx/60, dy=Ly/20, dz=Lz`：装置尺寸确定后自动派生；
- 第一阶段全部单元使用同一实测 `phi` 与 `k`，或同一实测常量各向异性张量；
- 第一阶段禁止随机渗透率/孔隙度场，地质非均质性后置；
- inlet/outlet manifold、压力测点与有效流动截面均由实验装置决定，不继承旧 point-well 几何；
- `nz=1` 仅代表厚度平均的 pseudo-3D slab，不能用于宣称解析真实三维垂向重力分异。

如果 slab 竖直且把模型 `y` 方向定义为物理竖直，可以研究 `ny=20` 网格上的**二维平面内浮力分异**；若重力方向垂直于 slab 平面，则 `nz=1` 无法解析厚度方向的饱和度/组分分层，必须改用 `nz>1` 的真正 3-D 模型。

旧 `benchmark_2d_common.hpp` 的 `1.20 m × 0.10 m × 0.10 m`、`phi=0.35`、`1500/150 mD` 仅保留为 legacy numerical benchmark，不属于新的实验设计。

**这些流动设置当前不执行。** 0D PVT 数值门禁已经通过，但最终进入流动还必须同时满足：PR binary calibration、CPA calibration/association、density validation、viscosity validation、laboratory porous-media/slab gate，以及 experiment-derived flow-control gate。统一依赖状态见 `pvt_acceptance/flow_entry_gate.csv`。

## 当前工作顺序

1. 保留已经通过的 `H2O_HEAVY_BINARY` production-kernel 数值筛查；优先追回/测量 360–380 °C、25–30 MPa 的 H2O–actual-Heavy 相数、相组成、密度和相界数据；
2. 在真实 Heavy calibration/hold-out 数据到位后，分别完成 PR 与 CPA 的 H2O–Heavy BIP/association 标定，不跨模型复制参数；随后再扩展到其它 H2O–lump；
3. 完成独立 density validation，并在另一数据层拟合必要的 volume translation；
4. 完成 intrinsic viscosity validation，再做 EOS-density + transport coupled validation；
5. 保持现已通过的 strict 0D PVT gate 作为所有后续参数更新的回归门禁；
6. 只有 `pvt_acceptance/flow_entry_gate.csv` 的所有硬依赖同时 PASS 后，才开始储层流动。
