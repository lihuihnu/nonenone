# 02 — 380 °C / 25 MPa 裂解产物原始数据集

## 目的

本阶段只建立可追溯的实验原始数据层，不进行 Light/Middle/Heavy 最终拟组分化，不生成 PR/CPA 参数，也不把跨样品数据、派生闭合值或图上坐标估读值伪装成主样品实验数据。

主样品继续使用 `01_RESEARCH_OBJECT.md` 锁定的 Zhao et al. (2023) 铜川 Chang 7 低成熟富有机质湖相页岩。

## 目标实验状态

- 温度：`380 °C` (`653.15 K`)
- 压力：`25 MPa`
- 反应时间：`4 h`
- 粒径：`75–106 μm`
- 两个 380 °C 重复实验：ESI experiments `5` and `6`
- Exp. 5：`21.99 g` shale + `22.01 g` water
- Exp. 6：`22.01 g` shale + `22.03 g` water

## RSC ESI 已追回：380 °C 精确重复实验数据

原论文 ESI Tables S1-S3 已解决此前“只能读图”的大部分缺口。原始重复值完整保存在：

`fluid_characterization/raw/rsc_esi_380c_replicates.csv`

### 油产率

Exp. 5：

- `52.27 mg/g shale`
- `346.05 mg/g TOC`

Exp. 6：

- `54.07 mg/g shale`
- `357.94 mg/g TOC`

两次重复的简单平均：`351.995 mg/g TOC`。

这与正文汇总值 `352.1 mg/g TOC` 一致到其报告精度。后续数据库同时保留正文值和重复实验原值，不能用平均值覆盖原始重复。

### 380 °C SARA

Exp. 5：Saturates `13.9%`、Aromatics `36.2%`、Resins `23.7%`、Asphaltenes `26.3%`。

Exp. 6：Saturates `12.5%`、Aromatics `38.2%`、Resins `24.5%`、Asphaltenes `24.8%`。

简单重复平均：Saturates `13.2%`、Aromatics `37.2%`、Resins `24.1%`、Asphaltenes `25.55%`。

这些区间是两个重复实验之间的 spread，不应被描述为仪器测量不确定度。

这一结果对后续模型非常关键：目标油在 380 °C SCW 转化后依然含有约一半的 `resins + asphaltenes`。因此最终 Heavy 绝不能在没有论证的情况下直接等同于一个非极性的单一正构烷烃或 squalane。

### 380 °C 产气

ESI Table S3 给出的原始体积产率为：

| component | Exp. 5 | Exp. 6 | unit |
|---|---:|---:|---|
| total gas | 5.80 | 6.21 | mL/g shale |
| H2 | 1.15 | 1.30 | mL/g shale |
| CO | 0 | 0 | mL/g shale |
| CH4 | 1.15 | 1.20 | mL/g shale |
| CO2 | 1.90 | 2.50 | mL/g shale |
| C2 | 0.41 | 0.59 | mL/g shale |
| C3 | 0.26 | 0.36 | mL/g shale |
| C4 | 0.12 | 0.20 | mL/g shale |
| C5 | 0.03 | 0.05 | mL/g shale |
| C6 | 0 | 0 | mL/g shale |

总气简单平均为 `6.005 mL/g shale`。以主样品 TOC `15.11 wt%` 进行派生归一化约为 `39.74 mL/g TOC`，但该 TOC 归一值是计算值而不是 ESI 原始字段。

两次实验给出的气体质量产率分别为 `5.58` 和 `7.34 mg/g shale`，简单平均 `6.46 mg/g shale`。

## 气体数据存在一个必须保留的统计口径问题

主论文正文报告 `H2 = 26.9%` at 380 °C，但直接把 ESI Table S3 的 `H2 mL/g shale` 除以同表的 `total gas mL/g shale`，两个 380 °C 重复实验都不能复现 `26.9%`。

因此本数据集采取以下规则：

1. `26.9%` 作为正文直接报告的 H2 composition statement 原样保存；
2. Table S3 的各气体 `mL/g shale` 作为原始 component-yield 数据原样保存；
3. 在弄清正文 Fig. 2c 的归一化/统计基准之前，不人为生成 CH4 / CO2 / C2 / C3+ 的百分比；
4. 后续若保留 Gas pseudo-component，优先从可追溯的 component yields 建立统一摩尔基准，而不是从图上读取百分比。

## 测量与回收边界

实验采用约 `80 mL` batch reactor（内径 `40 mm`、高度 `65 mm`），温度精度约 `±0.5 °C`，压力精度约 `±0.05 MPa`。

生成气使用 gas bag 收集、wet flowmeter 定量以及 Agilent 7890A GC（TCD + FID）分析。生成油用 CS2 从反应后页岩中洗脱，在约 `46 °C` 蒸发溶剂，并使用 CHNS elemental analysis 与 IATROSCAN MK-6 TLC 做 SARA。

原作者明确指出，46 °C 溶剂蒸发过程不可避免地损失沸点低于约 46 °C 的轻烃。因此：

> 回收液态油组成不是完整 C1+ 总产物组成。

后续若建立完整 Gas/Light/Middle/Heavy 初始流体，必须在同一质量/摩尔基准上联合气体与回收液体，并把低沸点损失作为实验不确定性和质量闭合问题处理。

## 另一个必须保持的物理边界

SARA 和气体数据来自 batch conversion 后的淬冷、降压和产品回收；它们不是 `380 °C / 25 MPa` 下 H2O + 裂解产物的原位平衡 tie-line。

因此：

- 本数据集用于定义预生成裂解产物的总体表征；
- PR/CPA 相平衡参数必须另用高温高压 H2O-hydrocarbon VLE/LLE/PVT 数据标定；
- 不允许直接用 post-quench SARA 去回归 380 °C / 25 MPa EOS 平衡组成。

## ACS 2023 配对纯干酪根原文转写包已审计

ACS DOI `10.1021/acs.iecr.3c02759` 研究的是由铜川 Chang 7 页岩酸洗并去除矿物后得到的 Type-II 干酪根。它与主样品具有很强的地质/材料配对意义，但处理状态不同，因此所有数据仍保持 `SECONDARY_PAIRED`，不能覆盖 RSC 完整 raw-shale 主样品数据。

本次上传的 Markdown 包由原始 11 页 PDF 转写，保留正文表格和原始 raster figures。Figure 6 自带数值标签，因此下列数值不是按坐标轴目测估计。

### 实验协议与材料分析

正文直接锁定：

- `25 MPa`；
- `300–700 °C` 温度系列；
- `2 h` 保温；
- 去离子水 : 干酪根质量比 `1:3`；
- `80 cm3` batch reactor；
- 温度精度 `±0.5 °C`；
- 压力精度 `±0.05 MPa`；
- 酸洗前页岩粒径 `120–180 μm`。

Table 2 给出的酸洗干酪根为：C `69.46 wt%`、H `6.01 wt%`、N `5.58 wt%`、S `2.39 wt%`、O `11.94 wt%`、moisture `2.71 wt%`、ash `1.91 wt%`、volatile matter `49.16 wt%`、fixed carbon `46.22 wt%`。完整表见 `fluid_characterization/raw/acs2023_table2_material_analysis.csv`。

### 380 °C 纯干酪根产油与 simulated distillation

正文报告 `380 °C / 25 MPa` 的纯干酪根产油峰值为 `0.19 g/g TOC`。

Figure 6 在 380 °C 的 simulated-distillation 四段质量分数为：

| fraction | boiling range | wt% |
|---|---|---:|
| Gasoline | IBP–180 °C | 0.81 |
| Diesel | 180–350 °C | 23.73 |
| Middle | 350–500 °C | 34.11 |
| Heavy | >500 °C | 41.35 |

四项合计 `100.00%`。完整的 300–500 °C 以及 free-oil Figure 6 数据已写入：

`fluid_characterization/raw/acs2023_figure6_distillation_sara.csv`

这组数据非常适合约束后续**候选** boiling-range lumping，但它仍是 pure-kerogen secondary paired data，不能直接当作完整 raw-shale 主样品最终 lump fractions。

### 380 °C SARA 的原图内部矛盾

Figure 6 在 380 °C 可直接读到：Saturates `7.21%`、Resins `44.45%`、Asphaltenes `16.78%`，但 Aromatics 的可见标签印为 `44.80%`。若全部照抄，合计为 `113.24%`，不可能是同一归一化 SARA 堆叠柱。

正文又明确给出该趋势起点 `saturates + aromatics ≈ 38.8%`，而 resins/asphaltenes 分别约 `44.5%` 和 `16.8%`。因此：

`Aromatics = 100 - 7.21 - 44.45 - 16.78 = 31.56%`

且 `7.21 + 31.56 = 38.77%`，与正文四舍五入的 `38.8%` 一致。

仓库中 `31.56%` 被明确标为 `DERIVED_MASS_CLOSURE_SOURCE_FIGURE_LABEL_INCONSISTENT`，绝不冒充原图直接标签。后续若取得 ACS SI 原始表，应以 SI 对这个源图问题作最终核对。

### 其它可直接使用的 ACS 原文数据

- Table 4，380 °C generated oil：`Xoxid=0.41`、`Xali=0.46`、`Xbrn=0.60`；完整 380–500 °C 系列见 `acs2023_table4_generated_oil_ftir_indices.csv`；
- Table 6，380 °C spent kerogen：C `79.58 wt%`、H `5.47 wt%`、N `5.06 wt%`、S `2.41 wt%`、H/C `0.82`；完整温度序列见 `acs2023_table6_spent_kerogen_ultimate.csv`；
- Section 3.2：free oil 的 C15–C25 为 `27.73%`；300–450 °C 生成油的 C15–C25 范围为 `17.36–39.51%`；500 °C 降为 `5.57%`；
- Section 3.1：650 °C gas yield `0.61 g/g TOC`，700 °C `0.88 g/g TOC`；
- Section 3.3：550 °C CO2 fraction `16.41%`，700 °C CO2 `20.62%`、CH4 `41.41%`；摘要给出 600 °C CH4 peak 约 `51%`、700 °C H2 约 `30%`。

摘要长表见 `fluid_characterization/raw/related_pure_kerogen_acs2023.csv`。

### ACS Supporting Information 仍未追回

原文说明 SI 包含 detailed SCW kerogen data、gas-production data、syngas-component data 和 generated-oil data。本次包不含完整 SI 数值表。因此 SI 仍有价值，主要用于：

- 核对 run-by-run product values；
- 解决 Figure 6 的 380 °C aromatic 标签矛盾；
- 获取正文图中未逐项标出的详细气体/产物数据。

但“ACS 380 °C 没有可用馏程数值”这一旧判断已经失效：**主文 Figure 6 已经给出了可直接录入的 exact visible labels**。

完整审计记录见 `fluid_characterization/raw/SUPPORTING_INFORMATION_AUDIT.md`。

## 同团队/同地层进一步搜索结果

### Xie et al. 2022, Oil Shale

DOI `10.3176/oil.2022.3.02` 使用 Ordos Basin `F317-181` 井样，油分析采用 Agilent 7890B chromatograph。论文指出生成油主要分布约 `C8-C56`，并在 380–450 °C 区间观察到低于 C16 的比例随温度升高下降、C16 以上比例增加。

这是非常有价值的候选 lump 边界证据，但它不是锁定的 Tongchuan outcrop 主样品：其样品来源和 TOC (`~16.25 wt%`) 均不同。因此它只能约束“怎样切”而不能决定“主样品各 lump 有多少”。

### 2023 Geoenergy Science and Engineering

DOI `10.1016/j.geoen.2023.211553` 属于同一研究方向，包含 360 °C/21 MPa 与 400 °C/25 MPa 等条件，可用于反应时间/温度趋势和未来对照工况设计。物理样品身份尚未与主样品逐项匹配，不并入主表。

### Lu et al. 2026

DOI `10.1016/j.jaap.2026.107757` 含 Chang-7 Type-II1 source-rock 样品，但其 380 °C 油产率为 `234.1 mg/g TOC`，与主样品 `~352 mg/g TOC` 差异明显。在拿到并逐项匹配样品表之前，不能认定为同一物理样品。

相关研究统一隔离存入 `fluid_characterization/raw/related_same_team_characterization.csv`。

## 当前文件

- `fluid_characterization/raw/README.md`
- `fluid_characterization/raw/primary_380c_25mpa_observations.csv`
- `fluid_characterization/raw/rsc_esi_380c_replicates.csv`
- `fluid_characterization/raw/related_pure_kerogen_acs2023.csv`
- `fluid_characterization/raw/acs2023_table2_material_analysis.csv`
- `fluid_characterization/raw/acs2023_figure6_distillation_sara.csv`
- `fluid_characterization/raw/acs2023_table4_generated_oil_ftir_indices.csv`
- `fluid_characterization/raw/acs2023_table6_spent_kerogen_ultimate.csv`
- `fluid_characterization/raw/SUPPORTING_INFORMATION_AUDIT.md`
- `fluid_characterization/raw/source_manifest.csv`
- `fluid_characterization/raw/data_gaps.csv`
- `fluid_characterization/raw/related_same_team_characterization.csv`

## 本阶段当前验收结论

### 已解决

- 380 °C 主样品两次重复实验的精确油产率、SARA、总产气、分组分体积产率和 mass balance；
- RSC ESI 的原始实验条件；
- ACS pure-kerogen 主文实验协议；
- ACS Table 2 / Table 4 / Table 6 的直接数值；
- ACS Figure 6 全温度 simulated-distillation/SARA 数值标签；
- ACS 380 °C paired pure-kerogen boiling-range profile；
- ACS Figure 6 380 °C SARA 源标签矛盾的显式隔离与闭合派生；
- 主样品与酸洗纯干酪根、同团队其他 Chang-7 / Ordos 数据之间的样品隔离规则。

### 仍未达到最终 pseudo-component 定值状态

剩余核心缺口按优先级为：

1. **同一完整 Tongchuan raw-shale 主样品 380 °C 生成油的 carbon-number distribution 或 simulated distillation**；
2. 回收油 density / specific gravity；
3. 回收油 average molecular weight；
4. full C1+ product mass/mole-basis reconstruction，包括低沸点液态烃回收损失；
5. RSC 正文 H2 26.9% 与 ESI Table S3 component yields 的归一化口径核对；
6. 高温高压油相黏度实验数据；
7. ACS 2023 detailed SI 数值表（仅作为 secondary paired constraint 与源图核对）。

因此下一阶段可以开始设计有文献约束的“候选 lumping 方案”，尤其可以把 ACS 380 °C 的四段馏程作为 secondary prior；但在没有同主样品馏程/碳数数据前，不应把候选切分和代表组分写成最终实验事实。
