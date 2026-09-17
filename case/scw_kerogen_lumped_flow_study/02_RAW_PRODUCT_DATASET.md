# 02 — 380 °C / 25 MPa 裂解产物原始数据集

## 目的

本阶段只建立可追溯的实验原始数据层，不进行 Light/Middle/Heavy 拟组分化，不生成 PR/CPA 参数，也不把图上目测值伪装成实验表格数据。

主样品继续使用 `01_RESEARCH_OBJECT.md` 锁定的 Zhao et al. (2023) 铜川 Chang 7 低成熟富有机质湖相页岩。

## 目标实验状态

- 温度：`380 °C` (`653.15 K`)
- 压力：`25 MPa`
- 反应时间：`4 h`
- 主温度序列水/页岩质量比：`1:1`
- 生烃/矿物效应实验粒径：`75–106 μm`

## 当前已取得的 A 级直接数据

### 产率

- 380 °C 油产率：`352.1 mg/g TOC`
- 这是该研究温度序列中的最高油产率。

380 °C 的精确总产气量在正文中没有逐项给数，论文明确将详细实验数据放在 ESI Table S2 / Fig. S1；在取得 ESI 数字前保持为空。

### 气体

主论文说明生成气主要包括：

- `H2`
- `CO2`
- `CH4`
- `C2`
- `C3+`

其中正文明确给出：

- `H2 = 26.9%` at 380 °C，且为 300–650 °C 序列中的局部最大值。

380 °C 下 `CH4 / C2 / C3+ / CO2` 的精确比例当前只在 Fig. 2c 图中存在，正文没有逐项数值，因此原始表保持空值并标记 `FIGURE_ONLY`。

### 液体油

论文对回收油进行了 SARA：

- Saturates
- Aromatics
- Resins
- Asphaltenes

Fig. 2b 给出温度序列，正文给出的可靠趋势是：350 °C 后随着温度继续升高，saturates 比例逐渐提高，而 asphaltenes / resins / aromatics 比例下降；但 380 °C 四项精确百分比当前未从机器可读表/ESI 取得，因此不得填入猜测值。

主论文未提供足以直接建立 `C6–C14 / C15–C20 / C21+` 的目标样品机器可读碳数分布或模拟蒸馏表。因此当前**不能**据此宣称 Light/Middle/Heavy 的最终切割范围与质量分数已经确定。

## 测量和回收边界

实验采用约 `80 mL` batch reactor（内径 `40 mm`、高度 `65 mm`），温度精度约 `±0.5 °C`，压力精度约 `±0.05 MPa`。

生成气：

- gas bag 收集；
- wet flowmeter 定量；
- Agilent 7890a GC，TCD + FID 分析。

生成油：

- CS2 从反应后页岩中洗脱；
- 约 `46 °C` 蒸发溶剂；
- CHNS elemental analysis；
- IATROSCAN MK-6 TLC 做 SARA。

原作者明确指出，46 °C 溶剂蒸发过程不可避免地损失沸点低于约 46 °C 的轻烃。因此：

> 回收液态油组成不是完整 C1+ 总产物组成。

后续若建立完整 Gas/Light/Middle/Heavy 初始流体，必须在同一质量基准上联合气体与液体信息，并把这部分轻烃损失作为实验不确定性处理。

## 另一个必须保持的物理边界

论文测得的 SARA 和气体组成是 batch conversion 完成后**淬冷、降压和回收后的产品分析**。

它不是：

`380 °C / 25 MPa 下 H2O + 裂解油的原位平衡相组成`。

因此以后：

- 本数据集用于定义“预生成裂解产物的总体组成/表征”；
- PR/CPA 相平衡参数必须另用高温高压 H2O–hydrocarbon VLE/LLE/PVT 数据标定；
- 不允许直接把 post-quench SARA 当成 380 °C / 25 MPa tie-line 去回归 EOS。

## 配对辅助数据

ACS 2023 `10.1021/acs.iecr.3c02759` 对酸洗去矿物的 Type-II kerogen 进行了 300–700 °C SCW 实验，并报告：

- 380 °C 纯干酪根油产率峰值 `0.19 g/g TOC`；
- 500 °C light distillates 约 `70%`；
- 500 °C asphaltene 约 `5%`；
- 600 °C CH4 比例峰值约 `51%`；
- 700 °C H2 约 `30%`；
- 700 °C gas yield `0.88 g/g TOC`；
- 对生成油使用了 SARA 和 simulated distillation。

这些数据被放入 `raw/related_pure_kerogen_acs2023.csv`，证据级别为配对辅助数据，不与主样品 raw-shale 380 °C 数据直接合并。

## 当前文件

- `fluid_characterization/raw/README.md`
- `fluid_characterization/raw/primary_380c_25mpa_observations.csv`
- `fluid_characterization/raw/source_manifest.csv`
- `fluid_characterization/raw/data_gaps.csv`
- `fluid_characterization/raw/related_pure_kerogen_acs2023.csv`

## 本阶段验收结论

### 已完成

- 锁定目标状态和实验协议；
- 建立带证据等级/数据状态的主样品原始观测表；
- 建立文献来源隔离规则；
- 建立关键缺口清单；
- 单独保存纯干酪根配对研究，而不混入主样品；
- 明确轻烃回收损失和 post-quench / in-situ composition 的物理区别。

### 尚未达到“可直接 lumping”状态

至少需要继续解决：

1. 380 °C 精确 SARA 四项数值；
2. 380 °C 精确总产气量；
3. 380 °C CH4 / CO2 / C2 / C3+ 精确比例；
4. 最关键：同一主样品 380 °C 生成油的碳数分布或模拟蒸馏数据；
5. 生成油平均分子量与 density/specific gravity；
6. 能够在统一质量基准上重建的 gas + recovered liquid composition。

在上述关键数据补齐之前，不进入最终 Light/Middle/Heavy 组成定值。