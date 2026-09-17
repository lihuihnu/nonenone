# 380 °C / 25 MPa Chang 7 SCW product raw dataset

本目录只保存能够追溯到公开实验来源的原始观测、明确标注的派生统计量、实验条件和数据缺口。它不是最终 EOS 输入表，也不是 Light/Middle/Heavy 拟组分表。

## 主数据对象

主样品固定为 Zhao et al. (2023, *Sustainable Energy & Fuels*, DOI `10.1039/D2SE01361D`) 报道的中国鄂尔多斯盆地南部长 7 段铜川露头低成熟富有机质湖相页岩。

主实验状态：

- `T = 380 °C = 653.15 K`
- `P = 25 MPa`
- `reaction time = 4 h`
- 380 °C ESI experiments `5` and `6`
- 水/页岩质量比约 `1:1`
- 页岩粒径 `75–106 μm`

## 数据状态

- `DIRECT_NUMERIC`: 正文、原始表格或补充材料明确给出的数值；
- `DIRECT_TABLE`: 原论文表格直接数值；
- `DIRECT_FIGURE_LABEL`: 原论文图中直接印刷的数值标签，不是人工按坐标估读；
- `DIRECT_QUALITATIVE`: 原文明确给出的定性信息；
- `DERIVED_FROM_ESI_DUPLICATES`: 仅由已保存的 ESI 重复实验作简单统计得到，不替代原始行；
- `DERIVED_MASS_CLOSURE_SOURCE_FIGURE_LABEL_INCONSISTENT`: 原图标签自相矛盾，仅在完整记录冲突后通过质量闭合得到的派生值；
- `DIRECT_NUMERIC_CONSISTENCY_CHECK_REQUIRED`: 原文直接报告，但与另一张原始表的统计口径尚需核对；
- `NOT_LOCATED_FOR_EXACT_SAMPLE`: 当前主论文和 ESI 审计没有找到同一物理样品对应的数据；
- `SECONDARY_PAIRED`: 同地质来源但处理状态不同的配对材料，只用于辅助约束；
- `CONTEXT_ONLY`: 不满足主样品身份条件，只能用于趋势或方法背景。

## 已经从 RSC ESI 补齐的 380 °C 主样品数据

- 两次重复油产率；
- 两次重复 SARA；
- 两次重复总产气量；
- H2 / CO / CH4 / CO2 / C2 / C3 / C4 / C5 / C6 的原始 `mL/g shale` 产率；
- 两次重复 mass balance；
- 两次重复样品/水装料量。

完整原始行见 `rsc_esi_380c_replicates.csv`。

## 已从上传的 ACS 2023 原文转写包补齐的 secondary paired 数据

论文：Zhao et al., *Industrial & Engineering Chemistry Research* 2023, DOI `10.1021/acs.iecr.3c02759`。

该研究使用的是**酸洗去矿物后的 Type-II 干酪根**，不是上面的完整 raw-shale 主样品，所以全部保持 `SECONDARY_PAIRED`。

当前已从原文表格、正文和带数值标签的 Figure 6 精确录入：

- Table 2：原始页岩与酸洗干酪根元素/工业分析；
- 实验协议：25 MPa、300–700 °C、2 h、水:干酪根质量比 1:3、80 cm3 反应釜；
- 380 °C 峰值产油：`0.19 g/g TOC`；
- Figure 6：300–500 °C 及 free oil 的 simulated-distillation 与 SARA 数值标签；
- 380 °C pure-kerogen simulated distillation：Gasoline `0.81%`、Diesel `23.73%`、Middle `34.11%`、Heavy `41.35%`；
- Table 4：380–500 °C 生成油 FT-IR indices；
- Table 6：初始/反应后干酪根元素组成与 H/C；
- 正文直接报告的 650/700 °C 产气和部分归一化气体组成锚点。

### Figure 6 的 380 °C SARA 源数据冲突

Figure 6 中 380 °C 的可见标签为 saturates `7.21%`、resins `44.45%`、asphaltenes `16.78%`，但 aromatic 标签印为 `44.80%`，四者会得到 `113.24%`，显然不闭合。正文同时给出该温度附近 `saturates + aromatics ≈ 38.8%`。

因此仓库：

- 保留这个原图内部矛盾；
- 不把 `44.80%` 当成可信的 380 °C aromatic 实验值；
- 仅以 `100 - 7.21 - 44.45 - 16.78 = 31.56%` 得到 aromatic 派生值，并标记为 `DERIVED_MASS_CLOSURE_SOURCE_FIGURE_LABEL_INCONSISTENT`；
- 未来若取得 ACS SI 原始表，以 SI 为优先核对来源。

## 仍然禁止的处理

1. 不允许把两个重复实验的平均值重新标成“实验原始值”。
2. 不允许因为正文报告 `H2 = 26.9%` 就用 RSC ESI component yields 强行归一出其余气体百分比；目前两套报告口径不能由简单相除完全复现。
3. 不允许把 ACS 酸洗纯干酪根实验、F317-181 井样、其他 Chang-7 样品直接并入主样品 380 °C / 25 MPa 数据行。
4. ACS pure-kerogen 的 380 °C simulated-distillation 可以用于设计候选 lump 边界/先验，但不能直接赋值给 intact-shale 主样品的最终 lump fractions。
5. 不允许将 SARA 直接转换成 C6–C14/C15–C20/C21+；需要碳数/馏程依据和明确映射规则。
6. 不允许用当前 `nC4/nC10/squalane` 筛选参数回填为实验原始数据。

## 文件

- `primary_380c_25mpa_observations.csv`: 主样品目标状态的摘要观测与明确标注的派生量；
- `rsc_esi_380c_replicates.csv`: RSC ESI 380 °C 两次重复实验的精确原始值；
- `related_pure_kerogen_acs2023.csv`: ACS 酸洗纯干酪根配对研究的摘要长表；
- `acs2023_table2_material_analysis.csv`: ACS Table 2 页岩/干酪根元素与工业分析；
- `acs2023_figure6_distillation_sara.csv`: ACS Figure 6 simulated-distillation 与 SARA 数值标签及一致性标记；
- `acs2023_table4_generated_oil_ftir_indices.csv`: ACS Table 4 生成油 FT-IR indices；
- `acs2023_table6_spent_kerogen_ultimate.csv`: ACS Table 6 初始/反应后干酪根元素分析；
- `SUPPORTING_INFORMATION_AUDIT.md`: RSC/ACS 原文、SI 与口径审计；
- `source_manifest.csv`: 文献来源、样品关系和允许用途；
- `data_gaps.csv`: 进入最终拟组分化前仍需补齐的关键缺口；
- `related_same_team_characterization.csv`: 同团队/同地层但不同或未确认同一物理样品的 GC、SARA 等证据。

## 当前结论

RSC ESI 已经把主样品 380 °C 的 SARA 和产气数据升级为可追溯的精确重复实验数据；本次 ACS 原文转写包又补出了一个很有价值的 **secondary paired 380 °C boiling-range constraint**，即 `0.81 / 23.73 / 34.11 / 41.35 wt%` 的 gasoline / diesel / middle / heavy 分布。

但最终 Light/Middle/Heavy 定值仍不能直接完成，核心缺口仍是：

1. **同一完整 Tongchuan raw-shale 主样品**的 carbon-number distribution 或 simulated distillation；
2. 回收油 density/specific gravity；
3. average molecular weight；
4. gas + recovered-liquid 的统一 C1+ 质量/摩尔基准及低沸点损失；
5. 高温高压黏度验证数据。
