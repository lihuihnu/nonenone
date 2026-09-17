# SCW–干酪根裂解产物拟组分流动算例集合

本目录用于后续建立“超临界水与预生成干酪根裂解产物的相平衡—黏度—多相多组分运移”数值实验集合。

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

因此当前烃相 lumping topology 直接采用这四个**实验馏程区间**。详细决策见 `03_EXPERIMENT_DRIVEN_LUMPING.md`，机器可读表见 `fluid_characterization/experimental_lumping_380c.csv`。

这组分数来自酸洗 Type-II 干酪根实验，属于 `SECONDARY_PAIRED`：它们是当前最直接的真实实验先验，但不能伪装成完整 raw-shale 主样品的同样分数。主样品如果后续获得同物理样品 simulated-distillation，应保留相同实验驱动原则并用主样品数据重新赋值。

## 当前流体拓扑

第一阶段非反应流动模型按以下结构组织：

- `H2O`
- `OIL_GASOLINE`：IBP–180 °C
- `OIL_DIESEL`：180–350 °C
- `OIL_MIDDLE`：350–500 °C
- `OIL_HEAVY`：>500 °C

气体产物 `H2 / CO2 / CH4 / C2+` 暂不强行并入上述油相 lumps。只有在 gas + recovered-liquid 的统一质量/摩尔基准建立、且不重复计算低沸点损失后，才增加独立气体组分或 gas lump。

## 代表组分和 EOS 参数全部后置

当前 **不指定任何 representative molecule**。尤其：

- `OIL_GASOLINE` 不等于 nC4；
- `OIL_DIESEL` 不等于 nC10；
- `OIL_HEAVY` 不等于 squalane。

每个实验馏程 lump 的 `MW`、density/SG、characterization boiling point、`Tc`、`Pc`、`omega`、PR/SW `kij`、CPA association parameters 和黏度参数，都必须在后续物性表征阶段基于实验数据或有文献依据的 petroleum-fraction correlation 单独确定。

## 当前数据审计状态

主样品仍以 Zhao et al. (*Sustainable Energy & Fuels*, 2023, DOI `10.1039/D2SE01361D`) 的完整 Chang 7 raw-shale `380 °C / 25 MPa / 4 h` 数据作为第一优先级，用于主样品油产率、SARA 和产气约束。

ACS 2023 酸洗纯干酪根数据用于提供直接的 boiling-range topology 和 paired prior。两套样品身份保持隔离，不能把 secondary paired fraction 重新标成 primary raw-shale measurement。

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

1. 为四个实验馏程 lump 补 `MW + SG/density + representative Tb`；
2. 生成并审计 petroleum-fraction `Tc/Pc/omega`；
3. 标定 H2O–lump 与 lump–lump PR/SW/CPA 相平衡参数；
4. 验证 380 °C / 25 MPa phase behavior；
5. 再进入黏度和流动算例。

在这些数据完成前，不用方便的纯组分替代物填空并宣称其为实验拟组分。
