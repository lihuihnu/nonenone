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

在任何 60×20×1 或其它储层流动算例开始前，必须分别完成：

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

目标温度为 `360 / 374 / 380 °C`，压力为 `25–30 MPa`。组成扫描以实验四个 oil lump 的摩尔比为中心，沿 H2O overall mole fraction 从 `0.01` 扫到 `0.995`，并加入 light-enriched / heavy-enriched 两个明确标为 deterministic sensitivity 的油组成族。

每个状态输出并检查：

- production flash convergence；
- final active-set stability；
- phase count / phase code；
- 两/三相组成；
- mass / molar density；
- LBC viscosity 与适用时的 IAPWS water viscosity reference；
- Oil / Gas / Water canonical phase role；
- material closure 与 active-phase fugacity closure。

主初始锚点固定为 `BASE oil ratio + z_H2O=0.20 + 25 MPa`，分别在 360/374/380 °C 下评估。PR 与 CPA 不要求预测相同相数；相数差异是模型结果。真正硬要求是**两个 EOS 各自都必须物理自洽**。

此外还输出 unrestricted O/G/W P–T map、oil/gas/water phase-onset envelope、restricted O/G bubble/dew projection，以及三温度下的 dense pressure-composition maps。

运行入口：

`make -C tools run-scw-kerogen-0d-pvt-acceptance`

硬门禁入口：

`make -C tools require-scw-kerogen-0d-pvt-acceptance`

即使 0D structural gate 通过，最终 flow comparison 仍必须同时满足 PR binary、CPA、density 和 viscosity 现有 gates；依赖关系见 `pvt_acceptance/flow_entry_gate.csv`。

当前严格 0D runtime gate 已经 **PASS**：

- PR registered scan：594/594；
- CPA registered scan：594/594；
- PR/CPA target-window P–T map health：全部 PASS；
- PR dense composition path：5103/5103；
- CPA dense composition path：5103/5103；
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

## 当前数据审计状态

主样品仍以 Zhao et al. (*Sustainable Energy & Fuels*, 2023, DOI `10.1039/D2SE01361D`) 的完整 Chang 7 raw-shale `380 °C / 25 MPa / 4 h` 数据作为第一优先级，用于主样品油产率、SARA 和产气约束。

ACS 2023 酸洗纯干酪根数据用于提供直接的 boiling-range topology、paired mass-fraction prior 和内部 carbon-number characterization。两套样品身份保持隔离，不能把 secondary paired 数据重新标成 primary raw-shale measurement。

## 储层数值实验（冻结）

计划中的伪三维单层结构仍保留为后续目标：`60 x 20 x 1`、均质岩石、等温全组分多相流、左侧注水右侧生产、380 °C / 25 MPa 主工况及 360 °C / 25 MPa 对照。

**这些流动设置当前不执行。** 0D PVT 数值门禁已经通过，但最终进入流动还必须同时满足：PR binary calibration、CPA calibration/association、density validation、viscosity validation，以及新的 laboratory porous-media gate。统一依赖状态见 `pvt_acceptance/flow_entry_gate.csv`。

## 当前工作顺序

1. 逐表追回并保存 H2O–lump 的高温高压 VLE/LLE/PVT 原始数据；
2. 分别完成 PR 与 CPA 的二元/association 标定，不跨模型复制参数；
3. 完成独立 density validation，并在另一数据层拟合必要的 volume translation；
4. 完成 intrinsic viscosity validation，再做 EOS-density + transport coupled validation；
5. 保持现已通过的 strict 0D PVT gate 作为所有后续参数更新的回归门禁；
6. 只有 `pvt_acceptance/flow_entry_gate.csv` 的所有硬依赖同时 PASS 后，才开始储层流动。
