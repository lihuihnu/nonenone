# SCW–干酪根裂解产物拟组分流动算例集合

本目录用于后续建立“超临界水与预生成干酪根裂解产物的相平衡—黏度—多相多组分运移”数值实验集合。

## 研究边界

- 主目标工况固定为 **380 °C、25 MPa**；后续保留 **360 °C、25 MPa** 作为同压力亚临界水对照候选。
- 研究对象是**预生成的干酪根裂解产物**与注入水之间的相平衡和流动耦合。
- 当前不把干酪根裂解反应、aquathermolysis、焦炭生成或能量方程纳入本算例集合。
- 拟组分边界、PR/SW/CPA 参数、二元作用参数和黏度模型必须有可追溯实验或文献依据，不以方便计算的代表分子反向定义真实流体。

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

第一阶段非反应流动模型按以下结构组织：

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

这些性质是**实验约束 + petroleum characterization correlation 的 provisional 值**，不是直接测得的临界性质。它们可以进入下一阶段 PR/SW/CPA screening，但最终参数必须通过高温高压相平衡/PVT 数据验证。

## Heavy 与 squalane 的边界

`OIL_HEAVY` 绝不默认等于 squalane。当前 >500 °C Heavy 的中心估计约为 `MW=660 g/mol`、`Tb=564 °C`、`SG=0.943`，并且实验 SARA 显示显著 resin/asphaltene character。Figure 5 可见 C38–C74 尾部也可能低估最不挥发/最极性的残余物，因此 Heavy 行明确标记为 `PROVISIONAL_HEAVY_TAIL_LOWER_BOUND_LIKE`。

squalane 仅保留为**重饱和烃 benchmark**，其参考数据独立保存在 `fluid_characterization/squalane_benchmark.csv`，不得把其纯组分性质复制到生产 `OIL_HEAVY` 行。

## 当前数据审计状态

主样品仍以 Zhao et al. (*Sustainable Energy & Fuels*, 2023, DOI `10.1039/D2SE01361D`) 的完整 Chang 7 raw-shale `380 °C / 25 MPa / 4 h` 数据作为第一优先级，用于主样品油产率、SARA 和产气约束。

ACS 2023 酸洗纯干酪根数据用于提供直接的 boiling-range topology、paired mass-fraction prior 和内部 carbon-number characterization。两套样品身份保持隔离，不能把 secondary paired 数据重新标成 primary raw-shale measurement。

## 计划的数值实验

基础几何保持伪三维单层结构：

- 网格：`60 x 20 x 1`；
- 均质岩石；
- 等温全组分多相流；
- 左侧注入水、右侧生产；
- 主工况：`380 °C / 25 MPa`；
- 后续控制：`360 °C / 25 MPa`；
- 重点比较温度跨临界变化，以及 PR、SW、CPA 的热力学结构差异。

主要观测量：各相组成/密度/相态、相黏度与相流度、四个实验馏程 lump 的生产组成和累计采出、轻重馏分选择性、质量守恒及数值收敛性。

## 后续工作顺序

1. 用当前 characterization 表建立 PR/SW 初始 pseudo-component 参数对象并做 380 °C / 25 MPa flash/stability preflight；
2. 搜集/标定 H2O–lump 与 lump–lump 高温高压 VLE/LLE/PVT 数据；
3. 单独审计 `OIL_HEAVY` 的 PNA/SARA 与 CPA association 表征；
4. 验证密度、volume shift 和黏度 closure；
5. 再进入 360/380 °C 对照流动算例。

整个过程中不允许用方便的纯组分替代物覆盖真实实验拟组分的 provenance。
