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
- `DIRECT_QUALITATIVE`: 原文明确给出的定性信息；
- `DERIVED_FROM_ESI_DUPLICATES`: 仅由已保存的 ESI 重复实验作简单统计得到，不替代原始行；
- `DIRECT_NUMERIC_CONSISTENCY_CHECK_REQUIRED`: 原文直接报告，但与另一张原始表的统计口径尚需核对；
- `NOT_LOCATED_FOR_EXACT_SAMPLE`: 当前主论文和 ESI 审计没有找到同一物理样品对应的数据；
- `SECONDARY_PAIRED`: 同团队相关材料或不同处理状态数据，只用于辅助约束；
- `CONTEXT_ONLY`: 不满足主样品身份条件，只能用于趋势或方法背景。

## 已经从 RSC ESI 补齐的 380 °C 数据

- 两次重复油产率；
- 两次重复 SARA；
- 两次重复总产气量；
- H2 / CO / CH4 / CO2 / C2 / C3 / C4 / C5 / C6 的原始 `mL/g shale` 产率；
- 两次重复 mass balance；
- 两次重复样品/水装料量。

完整原始行见 `rsc_esi_380c_replicates.csv`。

## 仍然禁止的处理

1. 不允许把两个重复实验的平均值重新标成“实验原始值”。
2. 不允许因为正文报告 `H2 = 26.9%` 就用 ESI component yields 强行归一出其余气体百分比；目前两套报告口径不能由简单相除完全复现。
3. 不允许把纯干酪根酸洗实验、F317-181 井样、其他 Chang-7 样品直接并入主样品 380 °C / 25 MPa 数据行。
4. 不允许将 SARA 直接转换成 C6–C14/C15–C20/C21+；需要同主样品蒸馏/碳数数据或明确的二级假设。
5. 不允许用当前 `nC4/nC10/squalane` 筛选参数回填为实验原始数据。

## 文件

- `primary_380c_25mpa_observations.csv`: 主样品目标状态的摘要观测与明确标注的派生量；
- `rsc_esi_380c_replicates.csv`: RSC ESI 380 °C 两次重复实验的精确原始值；
- `SUPPORTING_INFORMATION_AUDIT.md`: RSC/ACS SI 获取情况及数据口径审计；
- `source_manifest.csv`: 文献来源、样品关系和允许用途；
- `data_gaps.csv`: 进入最终拟组分化前仍需补齐的关键缺口；
- `related_pure_kerogen_acs2023.csv`: 酸洗纯干酪根 SCW 配对研究；
- `related_same_team_characterization.csv`: 同团队/同地层但不同或未确认同一物理样品的 GC、SARA 等证据。

## 当前结论

RSC ESI 已经把主样品 380 °C 的 SARA 和产气数据从“图中趋势”升级为可追溯的精确重复实验数据。当前真正阻塞最终 Light/Middle/Heavy 定值的核心问题已收缩为：

1. 同一 Tongchuan 主样品的 carbon-number distribution 或 simulated distillation；
2. 回收油 density/specific gravity；
3. average molecular weight；
4. gas + recovered-liquid 的统一 C1+ 质量/摩尔基准及低沸点损失；
5. 高温高压黏度验证数据。
