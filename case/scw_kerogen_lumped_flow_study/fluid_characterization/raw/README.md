# 380 °C / 25 MPa Chang 7 SCW product raw dataset

本目录只保存能够追溯到公开实验来源的原始观测、实验条件和数据缺口。它不是最终 EOS 输入表，也不是 Light/Middle/Heavy 拟组分表。

## 主数据对象

主样品固定为 Zhao et al. (2023, Sustainable Energy & Fuels, DOI 10.1039/D2SE01361D) 报道的中国鄂尔多斯盆地南部长 7 段铜川露头低成熟富有机质湖相页岩。

主实验状态：

- `T = 380 °C = 653.15 K`
- `P = 25 MPa`
- `reaction time = 4 h`
- 主温度序列水/页岩质量比 `1:1`
- 用于生烃特征分析的页岩粒径 `75–106 μm`

## 数据状态

每个字段使用以下状态之一：

- `DIRECT_NUMERIC`: 正文、表格或可访问补充材料明确给出的数值；
- `DIRECT_QUALITATIVE`: 原文明确给出的定性结论，但没有逐项数值；
- `FIGURE_ONLY`: 原论文图中存在，但当前未从原始表/ESI取得精确机器可读数值；
- `PENDING_PRIMARY_ESI`: 论文明确指出 ESI 有详细数据，但当前尚未获得对应原始数值；
- `SECONDARY_PAIRED`: 同一研究团队、相关材料/处理条件的数据，只用于辅助约束，不与主样品直接合并；
- `NOT_MEASURED_OR_NOT_REPORTED`: 当前来源没有测量或没有报告。

## 禁止事项

1. 不允许通过目测柱状图填入“精确实验值”而不标注数字化误差。
2. 不允许把纯干酪根酸洗实验、其他 Chang 7 样品或普通热解实验直接并入主样品 380 °C / 25 MPa 数据行。
3. 不允许将 SARA 直接转换成 C6–C14/C15–C20/C21+ 而没有额外蒸馏/GC 数据。
4. 不允许用当前 `nC4/nC10/squalane` 筛选参数回填为实验原始数据。

## 文件

- `primary_380c_25mpa_observations.csv`: 主样品在目标状态下当前可核实的原始观测；
- `source_manifest.csv`: 文献来源与用途；
- `data_gaps.csv`: 进入拟组分化之前必须继续补齐的关键缺口；
- `related_pure_kerogen_acs2023.csv`: 同团队纯干酪根 SCW 实验，仅作为配对辅助数据，不与主样品直接混合。

## 当前结论

目前能直接锁定的是实验状态、油产率、部分气体信息、测量方法和轻组分损失边界；完整 380 °C SARA 数字、380 °C 总产气量、CH4/C2/C3+ 精确比例，以及真正用于碳数 lumping 的模拟蒸馏/GC 分布仍需继续获取原始补充数据或经过受控数字化后才能进入下一阶段。