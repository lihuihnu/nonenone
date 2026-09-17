# 02 — 380 °C / 25 MPa 裂解产物原始数据集

## 目的

本阶段只建立可追溯的实验原始数据层，不进行 Light/Middle/Heavy 最终拟组分化，不生成 PR/CPA 参数，也不把跨样品数据或图上目测值伪装成主样品实验数据。

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

两次重复的简单平均：

- `351.995 mg/g TOC`

这与正文汇总值 `352.1 mg/g TOC` 一致到其报告精度。后续数据库同时保留正文值和重复实验原值，不能用平均值覆盖原始重复。

### 380 °C SARA

Exp. 5：

- Saturates `13.9%`
- Aromatics `36.2%`
- Resins `23.7%`
- Asphaltenes `26.3%`

Exp. 6：

- Saturates `12.5%`
- Aromatics `38.2%`
- Resins `24.5%`
- Asphaltenes `24.8%`

简单重复平均：

- Saturates `13.2%`
- Aromatics `37.2%`
- Resins `24.1%`
- Asphaltenes `25.55%`

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

主论文正文报告：

- `H2 = 26.9%` at 380 °C.

但直接把 ESI Table S3 的 `H2 mL/g shale` 除以同表的 `total gas mL/g shale`，两个 380 °C 重复实验都不能复现 `26.9%`。

因此本数据集采取以下规则：

1. `26.9%` 作为正文直接报告的 H2 composition statement 原样保存；
2. Table S3 的各气体 `mL/g shale` 作为原始 component-yield 数据原样保存；
3. 在弄清正文 Fig. 2c 的归一化/统计基准之前，不人为生成 CH4 / CO2 / C2 / C3+ 的百分比；
4. 后续若保留 Gas pseudo-component，优先从可追溯的 component yields 建立统一摩尔基准，而不是从图上读取百分比。

## 测量与回收边界

实验采用约 `80 mL` batch reactor（内径 `40 mm`、高度 `65 mm`），温度精度约 `±0.5 °C`，压力精度约 `±0.05 MPa`。

生成气：

- gas bag 收集；
- wet flowmeter 定量；
- Agilent 7890A GC，TCD + FID 分析。

生成油：

- CS2 从反应后页岩中洗脱；
- 约 `46 °C` 蒸发溶剂；
- CHNS elemental analysis；
- IATROSCAN MK-6 TLC 做 SARA。

原作者明确指出，46 °C 溶剂蒸发过程不可避免地损失沸点低于约 46 °C 的轻烃。因此：

> 回收液态油组成不是完整 C1+ 总产物组成。

后续若建立完整 Gas/Light/Middle/Heavy 初始流体，必须在同一质量/摩尔基准上联合气体与回收液体，并把低沸点损失作为实验不确定性和质量闭合问题处理。

## 另一个必须保持的物理边界

SARA 和气体数据来自 batch conversion 后的淬冷、降压和产品回收；它们不是：

`380 °C / 25 MPa 下 H2O + 裂解产物的原位平衡 tie-line`。

因此：

- 本数据集用于定义预生成裂解产物的总体表征；
- PR/CPA 相平衡参数必须另用高温高压 H2O-hydrocarbon VLE/LLE/PVT 数据标定；
- 不允许直接用 post-quench SARA 去回归 380 °C / 25 MPa EOS 平衡组成。

## ACS 2023 配对纯干酪根正文已审计

ACS DOI `10.1021/acs.iecr.3c02759` 研究的是由铜川 Chang 7 页岩酸洗并去除矿物后得到的 Type-II 干酪根。它与主样品具有很强的地质/材料配对意义，但处理状态不同，因此所有数据仍隔离保存在：

`fluid_characterization/raw/related_pure_kerogen_acs2023.csv`

正文直接锁定的实验协议包括：

- `25 MPa`；
- `300–700 °C` 温度系列；
- `2 h` 保温；
- 去离子水 : 干酪根质量比 `1:3`；
- `80 cm3` batch reactor；
- 温度精度 `±0.5 °C`；
- 压力精度 `±0.05 MPa`；
- 酸洗前页岩粒径 `120–180 μm`。

与当前 380 °C 流体表征最相关的正文数值是：

- `380 °C / 25 MPa`：纯干酪根产油峰值 `0.19 g/g TOC`；
- `500 °C / 25 MPa`：light distillates `70%`，asphaltene `5%`；
- `600 °C / 25 MPa`：CH4 fraction peak `51%`；
- `700 °C / 25 MPa`：H2 fraction `30%`，gas yield `0.88 g/g TOC`。

正文还明确使用 SARA 与 simulated distillation 表征生成油，并指出汽油/柴油等轻质馏分随温度提高而增加。这些是有用的趋势约束，但不能反推出主样品 380 °C 的 Light/Middle/Heavy 比例。

### ACS Supporting Information 仍未追回

出版页面确认 SI 包含 detailed SCW kerogen data、gas-production data、syngas-component data 和 generated-oil data。当前已经完成**主文**审计，但尚未取得完整 SI 数值表，因此：

- 不填任何未直接报告的 380 °C simulated-distillation 百分比；
- 不从图上目测并伪装成精确值；
- 后续取得 SI 后单独录入，并继续保持 `SECONDARY_PAIRED` 身份；
- 即使取得，它也不能替代完整 raw-shale 主样品数据。

完整审计记录见：

`fluid_characterization/raw/SUPPORTING_INFORMATION_AUDIT.md`

## 同团队/同地层进一步搜索结果

### Xie et al. 2022, Oil Shale

DOI `10.3176/oil.2022.3.02` 使用 Ordos Basin `F317-181` 井样，油分析采用 Agilent 7890B chromatograph。论文指出生成油主要分布约 `C8-C56`，并在 380–450 °C 区间观察到低于 C16 的比例随温度升高下降、C16 以上比例增加。

这是非常有价值的候选 lump 边界证据，但它不是锁定的 Tongchuan outcrop 主样品：其样品来源和 TOC (`~16.25 wt%`) 均不同。因此它只能约束“怎样切”而不能决定“主样品各 lump 有多少”。

### 2023 Geoenergy Science and Engineering

DOI `10.1016/j.geoen.2023.211553` 属于同一研究方向，包含 360 °C/21 MPa 与 400 °C/25 MPa 等条件，可用于反应时间/温度趋势和未来对照工况设计。物理样品身份尚未与主样品逐项匹配，不并入主表。

### Lu et al. 2026

DOI `10.1016/j.jaap.2026.107757` 含 Chang-7 Type-II1 source-rock 样品，但其 380 °C 油产率为 `234.1 mg/g TOC`，与主样品 `~352 mg/g TOC` 差异明显。在拿到并逐项匹配样品表之前，不能认定为同一物理样品。

相关研究统一隔离存入：

`fluid_characterization/raw/related_same_team_characterization.csv`

## 当前文件

- `fluid_characterization/raw/README.md`
- `fluid_characterization/raw/primary_380c_25mpa_observations.csv`
- `fluid_characterization/raw/rsc_esi_380c_replicates.csv`
- `fluid_characterization/raw/SUPPORTING_INFORMATION_AUDIT.md`
- `fluid_characterization/raw/source_manifest.csv`
- `fluid_characterization/raw/data_gaps.csv`
- `fluid_characterization/raw/related_pure_kerogen_acs2023.csv`
- `fluid_characterization/raw/related_same_team_characterization.csv`

## 本阶段当前验收结论

### 已解决

- 380 °C 两次重复实验的精确油产率；
- 380 °C 精确 SARA；
- 380 °C 精确总产气量；
- 380 °C H2 / CH4 / CO2 / C2 / C3 / C4 / C5 / C6 原始体积产率；
- 380 °C 两次重复实验 mass balance (`90%`, `91%`)；
- RSC ESI 的原始实验条件；
- ACS 纯干酪根研究的主文实验协议与直接数值锚点；
- 主样品与酸洗纯干酪根、同团队其他 Chang-7 / Ordos 数据之间的样品隔离规则。

### 仍未达到最终 pseudo-component 定值状态

剩余核心缺口按优先级为：

1. **同一 Tongchuan 主样品 380 °C 生成油的 carbon-number distribution 或 simulated distillation**；
2. 回收油 density / specific gravity；
3. 回收油 average molecular weight；
4. full C1+ product mass/mole-basis reconstruction，包括低沸点液态烃回收损失；
5. RSC 正文 H2 26.9% 与 ESI Table S3 component yields 的归一化口径核对；
6. 高温高压油相黏度实验数据；
7. ACS 2023 detailed SI 数值表（仅作为 secondary paired constraint）。

因此下一阶段可以开始设计“候选 lumping 方案”，但在没有同主样品碳数/馏程数据前，不应把候选切分和代表组分写成最终实验事实。
