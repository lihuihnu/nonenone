# 02 — 380 °C / 25 MPa 裂解产物原始数据集

> 本文件只负责原始/派生数据审计。正式 lumping 决策已经迁移到 `03_EXPERIMENT_DRIVEN_LUMPING.md`：当前 hydrocarbon topology 直接采用实验 Figure 6 的四段馏程 `IBP–180 / 180–350 / 350–500 / >500 °C`，不再由 `nC4/nC10/squalane` 或预设碳数区间定义。机器可读定义见 `fluid_characterization/experimental_lumping_380c.csv`。

## 目的

本阶段只建立可追溯的实验原始数据层，不生成 PR/CPA 参数，也不把跨样品数据、派生闭合值或图上坐标估读值伪装成主样品实验数据。

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

原始重复值保存在 `fluid_characterization/raw/rsc_esi_380c_replicates.csv`。

### 油产率

Exp. 5：`52.27 mg/g shale`，`346.05 mg/g TOC`。

Exp. 6：`54.07 mg/g shale`，`357.94 mg/g TOC`。

两次重复简单平均为 `351.995 mg/g TOC`，与正文汇总 `352.1 mg/g TOC` 一致到报告精度。平均值只作为派生统计量，不覆盖原始重复。

### 380 °C SARA

Exp. 5：Saturates `13.9%`、Aromatics `36.2%`、Resins `23.7%`、Asphaltenes `26.3%`。

Exp. 6：Saturates `12.5%`、Aromatics `38.2%`、Resins `24.5%`、Asphaltenes `24.8%`。

简单平均：Saturates `13.2%`、Aromatics `37.2%`、Resins `24.1%`、Asphaltenes `25.55%`。

该结果说明重质部分具有显著 resin/asphaltene 特征，因此 `OIL_HEAVY` 不能在没有独立验证时等同于 squalane 或单一非极性正构烷烃。

### 380 °C 产气

ESI Table S3 原始体积产率：

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

总气简单平均 `6.005 mL/g shale`。按主样品 TOC `15.11 wt%` 派生归一约 `39.74 mL/g TOC`，但这不是 ESI 原始字段。气体质量产率为 `5.58` 和 `7.34 mg/g shale`。

## 气体统计口径问题

主论文正文报告 `H2 = 26.9%` at 380 °C，但直接用 ESI Table S3 的 H2 yield / total-gas yield 不能复现该比例。因此：

1. `26.9%` 原样保存为正文 composition statement；
2. Table S3 component yields 原样保存；
3. 在归一化基准澄清前不制造其它气体百分比；
4. 气体组分不得自动折入四个 recovered-oil boiling-range lumps。

## 测量与回收边界

RSC 实验采用约 `80 mL` batch reactor，温度精度约 `±0.5 °C`，压力精度约 `±0.05 MPa`。生成气采用 gas bag、wet flowmeter 和 Agilent 7890A GC（TCD + FID）；生成油用 CS2 洗脱，并在约 `46 °C` 蒸发溶剂后做 CHNS 和 SARA。

原作者明确指出低于约 46 °C 沸点的轻烃会因溶剂蒸发而损失，因此回收液态油不是完整 C1+ 总产物。gas + recovered liquid 必须在共同质量/摩尔基准上重构后，才能建立含气体产物的完整 compositional feed。

SARA 和气体数据也不是 `380 °C / 25 MPa` 原位 equilibrium tie-line，因此不得直接用于回归 EOS 平衡组成。

## ACS 2023 配对纯干酪根原文数据

ACS DOI `10.1021/acs.iecr.3c02759` 使用由铜川 Chang 7 页岩酸洗去矿物后得到的 Type-II kerogen。该数据保持 `SECONDARY_PAIRED`，但它提供了目前最直接的真实 380 °C generated-oil distillation evidence。

实验协议：`25 MPa`、`300–700 °C`、`2 h`、水:干酪根质量比 `1:3`、`80 cm3` reactor、温度精度 `±0.5 °C`、压力精度 `±0.05 MPa`。

Table 2、Table 4、Table 6 的完整直接数据分别保存在：

- `fluid_characterization/raw/acs2023_table2_material_analysis.csv`
- `fluid_characterization/raw/acs2023_table4_generated_oil_ftir_indices.csv`
- `fluid_characterization/raw/acs2023_table6_spent_kerogen_ultimate.csv`

### 380 °C simulated distillation：用于正式确定 lump topology

Figure 6 直接标签：

| lump | boiling range | wt% of recovered oil |
|---|---|---:|
| `OIL_GASOLINE` | IBP–180 °C | 0.81 |
| `OIL_DIESEL` | 180–350 °C | 23.73 |
| `OIL_MIDDLE` | 350–500 °C | 34.11 |
| `OIL_HEAVY` | >500 °C | 41.35 |

四项合计 `100.00%`。这些实验切分现已定义本算例的 hydrocarbon lump boundaries；不再把人为碳数范围作为并列候选。

完整 300–500 °C 及 free-oil Figure 6 数据见 `fluid_characterization/raw/acs2023_figure6_distillation_sara.csv`。

需要保持样品层级：上述 380 °C fraction weights 是 pure-kerogen secondary paired measurements；若后续获得同一 intact raw-shale 主样品的 simulated distillation，则用主样品在**同四个实验馏程区间**中的权重替换当前 paired prior，而不是重新按方便的代表分子切 lump。

### 380 °C SARA 原图矛盾

Figure 6 可见标签为 Saturates `7.21%`、Resins `44.45%`、Asphaltenes `16.78%`，Aromatics 印为 `44.80%`，总和为 `113.24%`，不闭合。正文同时给出 saturates + aromatics ≈ `38.8%`。

因此仓库只把 `Aromatics = 31.56%` 保存为显式质量闭合派生值：`100 - 7.21 - 44.45 - 16.78 = 31.56%`，状态为 `DERIVED_MASS_CLOSURE_SOURCE_FIGURE_LABEL_INCONSISTENT`。未来若取得 ACS SI 原始表，用 SI 核对该源图问题。

### 其它 ACS 原文锚点

- 380 °C pure-kerogen oil yield：`0.19 g/g TOC`；
- Table 4 at 380 °C：`Xoxid=0.41`、`Xali=0.46`、`Xbrn=0.60`；
- Table 6 spent kerogen at 380 °C：C `79.58 wt%`、H `5.47 wt%`、N `5.06 wt%`、S `2.41 wt%`、H/C `0.82`；
- free oil C15–C25：`27.73%`；generated oil C15–C25 across 300–450 °C：`17.36–39.51%`；500 °C：`5.57%`；
- gas yield：650 °C `0.61 g/g TOC`、700 °C `0.88 g/g TOC`。

## 当前结论

### 已经确定

- hydrocarbon lumping topology：`IBP–180 / 180–350 / 350–500 / >500 °C`；
- paired pure-kerogen 380 °C recovered-oil weights：`0.81 / 23.73 / 34.11 / 41.35 wt%`；
- 不再以 `nC4/nC10/squalane` 或任意碳数分段定义真实流体。

### 仍待表征

1. intact raw-shale 主样品在上述四个馏程区间中的准确权重；
2. 每个 lump 的 MW、density/SG 和 characterization Tb；
3. `Tc/Pc/omega` 与 PR/SW/CPA 参数；
4. gas + liquid 统一 C1+ feed reconstruction；
5. 380 °C / 25 MPa 高温高压 VLE/LLE/PVT 与黏度验证。

因此下一步不是再次选择 lump 边界，而是**对已经由实验确定的四个馏程 lump 做物性表征**。
