# 380 °C / 25 MPa Chang 7 SCW product raw dataset

本目录只保存能够追溯到公开实验来源的原始观测、明确标注的派生统计量、实验条件和数据缺口。它不是最终 EOS 输入表。

## 主数据对象

主样品固定为 Zhao et al. (2023, *Sustainable Energy & Fuels*, DOI `10.1039/D2SE01361D`) 报道的中国鄂尔多斯盆地南部长 7 段铜川露头低成熟富有机质湖相页岩。

主实验状态：`380 °C / 25 MPa / 4 h`，380 °C ESI experiments 5 and 6，水/页岩质量比约 1:1，页岩粒径 75–106 μm。

## 数据状态

- `DIRECT_NUMERIC` / `DIRECT_TABLE`: 正文、表格或补充材料直接数值；
- `DIRECT_FIGURE_LABEL`: 原论文图中直接印刷的数值标签，不是按坐标估读；
- `DIRECT_QUALITATIVE`: 原文定性信息；
- `DERIVED_FROM_ESI_DUPLICATES`: 由已保存重复实验作简单统计得到；
- `DERIVED_MASS_CLOSURE_SOURCE_FIGURE_LABEL_INCONSISTENT`: 源图自相矛盾时才使用的显式派生闭合值；
- `SECONDARY_PAIRED`: 同地质来源但处理状态不同的配对材料；
- `CONTEXT_ONLY`: 只允许用于趋势或方法背景。

## RSC 主样品数据

RSC ESI 已提供 380 °C 两次重复的油产率、SARA、总产气量、H2/CO/CH4/CO2/C2-C6 体积产率、mass balance 以及装料量。原始行见 `rsc_esi_380c_replicates.csv`。

## ACS 2023 配对纯干酪根数据

Zhao et al., *Industrial & Engineering Chemistry Research* 2023, DOI `10.1021/acs.iecr.3c02759` 使用酸洗去矿物后的 Type-II 干酪根，因此保持 `SECONDARY_PAIRED`。

已从原文直接录入：Table 2 材料分析、实验协议、380 °C 峰值产油、Figure 6 simulated distillation/SARA、Table 4 FT-IR indices、Table 6 spent-kerogen ultimate analysis，以及正文直接报告的部分气体数据。

### 380 °C Figure 6 直接决定当前 lumping topology

380 °C recovered generated oil 的 simulated-distillation 直接标签为：

- Gasoline, IBP–180 °C: `0.81 wt%`
- Diesel, 180–350 °C: `23.73 wt%`
- Middle distillate, 350–500 °C: `34.11 wt%`
- Heavy residue, >500 °C: `41.35 wt%`

四项闭合到 `100.00 wt%`。

因此本项目后续的 pseudo-component **边界**不再由 `nC4/nC10/squalane` 或人为碳数区间预先规定，而直接采用上述真实实验馏程切分。正式决策见 `../../03_EXPERIMENT_DRIVEN_LUMPING.md`；机器可读表见 `../experimental_lumping_380c.csv`。

注意：这四个比例是 pure-kerogen secondary-paired measurements。它们可以作为当前 paired experimental prior，但不能重新标成 intact raw-shale primary fractions。若未来取得主样品同物理样品 simulated distillation，应在相同实验驱动原则下更新四个 fraction weights，而不是重新回到方便计算的代表分子切分。

### Figure 6 的 380 °C SARA 源数据冲突

Figure 6 中 380 °C 的可见标签为 saturates `7.21%`、resins `44.45%`、asphaltenes `16.78%`，但 aromatic 标签印为 `44.80%`，总和变成 `113.24%`。正文同时给出 `saturates + aromatics ≈ 38.8%`。

仓库保留源图冲突，并仅以质量闭合得到 aromatics `31.56%`，标记为 `DERIVED_MASS_CLOSURE_SOURCE_FIGURE_LABEL_INCONSISTENT`。若以后取得 ACS SI 原始表，以 SI 核对。

## 仍然禁止的处理

1. 不把重复实验平均值重新标成原始实验值。
2. 不用 RSC 的 H2 百分比强行归一制造其它气体百分比。
3. 不把 ACS pure-kerogen、F317-181 或其它 Chang-7 样品直接并入 primary raw-shale 行。
4. 不从 SARA 反推 carbon-number lumps。
5. 不用 `nC4/nC10/squalane` 或其它方便的纯组分反向定义实验 lump。
6. 不在缺失 MW/SG/Tb/Tc/Pc/omega 时用代表正构烷烃填空并称为实验物性。

## 关键文件

- `primary_380c_25mpa_observations.csv`
- `rsc_esi_380c_replicates.csv`
- `related_pure_kerogen_acs2023.csv`
- `acs2023_table2_material_analysis.csv`
- `acs2023_figure6_distillation_sara.csv`
- `acs2023_table4_generated_oil_ftir_indices.csv`
- `acs2023_table6_spent_kerogen_ultimate.csv`
- `SUPPORTING_INFORMATION_AUDIT.md`
- `source_manifest.csv`
- `data_gaps.csv`
- `../experimental_lumping_380c.csv`

## 当前结论

**Lumping topology 已经可以由真实实验确定**：IBP–180 / 180–350 / 350–500 / >500 °C 四段。

现在真正阻塞后续 EOS/transport 的不是“该怎么切 lump”，而是每个实验馏程 fraction 的 `MW`、density/SG、characterization Tb、`Tc/Pc/omega`、CPA association 信息、高温高压相平衡和黏度数据，以及完整 raw-shale 主样品在这四段中的准确 fraction weights。
