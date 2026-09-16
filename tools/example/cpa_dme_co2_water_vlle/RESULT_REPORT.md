# CO2–DME–H2O CPA 三相闪蒸结果

## 结论

Our CPA 在 308.15 K、19.0–58.8 bar 的 4 个公开三元 VLLE 状态点上均
自动返回 DME-rich liquid、vapor 和 water-rich liquid 三个相。所有点的
相态码均为 `7`，总体组成重构误差不超过 `8.89e-15`，三相逸度对数最大差
不超过 `1.24e-10`。因此，这组计算通过了三相闪蒸的相态识别、物料守恒和
相平衡数值条件。

与实验的 36 个相组成值相比，总体 MAE 为 `0.01387`（摩尔分数），最大绝对
误差为 `0.04156`。主要组分和 DME、CO2 的预测与 Folas 已发表的 CPA 误差
基本一致；痕量水的相对误差很大，但其绝对误差仍小。这说明该算例同时复现了
公开 CPA 模型的优势与已知局限，不能把痕量组分 AARD 误读为闪蒸数值失败。

## 实验与计算设置

实验数据来自 Laursen、Rasmussen 和 Andersen 的 CO2–DME–H2O VLLE
测量。温度为 308.15 K，压力分别为 19.0、31.7、46.0 和 58.8 bar；每个点
都报告 water-rich liquid、DME-rich liquid 和 vapor 的完整组成。

计算使用 Folas 给出的 simplified CPA、DME 纯组分参数、4C water 参数、
三个二元相互作用参数和 water–DME mCR-1 交叉缔合参数。为避免用结果反拟合，
所有参数在运行前固定，未按本次误差调整。

原论文未报告总体进料。每个点的总体组成取三组实验相组成的等权平均：

| P (bar) | z(H2O) | z(DME) | z(CO2) |
|---:|---:|---:|---:|
| 19.0 | 0.32787 | 0.38827 | 0.28387 |
| 31.7 | 0.32460 | 0.25850 | 0.41690 |
| 46.0 | 0.33337 | 0.15553 | 0.51110 |
| 58.8 | 0.33957 | 0.08580 | 0.57463 |

该总体组成位于公开三相组成形成的三角形内部，保证三个相都有正相分率。相分率
是这个人工进料的结果，不与论文对比；真正用于实验对比的是三个平衡相的组成。

## 三相数值结果

`beta_oil` 对应 DME-rich liquid，`beta_gas` 对应 vapor，`beta_water`
对应 water-rich liquid。

| P (bar) | 迭代 | beta_oil | beta_gas | beta_water | 物料闭合 | 最大 log-fugacity 差 |
|---:|---:|---:|---:|---:|---:|---:|
| 19.0 | 14 | 0.29835 | 0.34266 | 0.35898 | 5.55e-17 | 1.24e-10 |
| 31.7 | 14 | 0.29690 | 0.35274 | 0.35036 | 5.55e-17 | 2.13e-14 |
| 46.0 | 13 | 0.28834 | 0.35784 | 0.35382 | 5.55e-17 | 3.05e-11 |
| 58.8 | 13 | 0.28244 | 0.36241 | 0.35514 | 8.88e-15 | 1.60e-12 |

压力升高时，实验与计算都表现出相同趋势：CO2 在 DME-rich liquid 和 vapor
中的摩尔分数增加，DME 相应减少；water-rich liquid 始终以 H2O 为主。三相
组成随压力的演化方向被正确恢复。

## 相组成对比

各单元格按“实验 / Our”列出摩尔分数：

| P (bar) | 相 | x(H2O) | x(DME) | x(CO2) |
|---:|---|---:|---:|---:|
| 19.0 | water-rich | 0.9058 / 0.8768 | 0.0862 / 0.1119 | 0.0080 / 0.0113 |
| 19.0 | DME-rich | 0.0686 / 0.0402 | 0.7142 / 0.7469 | 0.2172 / 0.2129 |
| 19.0 | vapor | 0.0092 / 0.0032 | 0.3644 / 0.3655 | 0.6264 / 0.6312 |
| 31.7 | water-rich | 0.9290 / 0.9085 | 0.0546 / 0.0737 | 0.0164 / 0.0178 |
| 31.7 | DME-rich | 0.0339 / 0.0187 | 0.5227 / 0.5490 | 0.4434 / 0.4323 |
| 31.7 | vapor | 0.0109 / 0.0021 | 0.1982 / 0.1976 | 0.7909 / 0.8004 |
| 46.0 | water-rich | 0.9549 / 0.9347 | 0.0288 / 0.0437 | 0.0163 / 0.0216 |
| 46.0 | DME-rich | 0.0350 / 0.0074 | 0.3207 / 0.3440 | 0.6443 / 0.6487 |
| 46.0 | vapor | 0.0102 / 0.0015 | 0.1171 / 0.1142 | 0.8727 / 0.8843 |
| 58.8 | water-rich | 0.9612 / 0.9526 | 0.0154 / 0.0241 | 0.0234 / 0.0234 |
| 58.8 | DME-rich | 0.0445 / 0.0029 | 0.1688 / 0.1852 | 0.7867 / 0.8119 |
| 58.8 | vapor | 0.0130 / 0.0012 | 0.0732 / 0.0688 | 0.9138 / 0.9300 |

## 图像解读

中英文组成 parity 图生成在 `results/figures/zh/` 和 `results/figures/en/`。
横轴为实验摩尔分数，纵轴为 Our CPA 计算摩尔分数；点越接近 1:1 虚线，实验与
计算越吻合。颜色与实心/空心/斜线填充共同区分三个相，圆形、方形、三角形分别
表示 H2O、DME、CO2。36 个点整体沿 1:1 线分布，偏离较明显的点主要位于低摩尔
分数区，对应痕量水和少量溶解组分；完整数值可由上一节表格读取。

## 与公开 CPA 误差对照

下表中的 Our AARD 只统计 4 个三个组分均非零的三元点。Folas 表 7.5 的
published AAD 使用其完整 308.15 K 数据集合，其中还包含 7.5 bar 的二元端点，
所以两列不是逐点完全相同的统计样本；它们用于检查量级和分相规律是否复现。

| 相 | 组分 | Our AARD (%) | Folas published AAD (%) |
|---|---|---:|---:|
| water-rich | H2O | 2.10 | 1.8 |
| water-rich | DME | 43.23 | 41.0 |
| water-rich | CO2 | 20.60 | 15.4 |
| DME-rich | H2O | 64.63 | 65.2 |
| DME-rich | DME | 6.64 | 6.6 |
| DME-rich | CO2 | 2.09 | 2.1 |
| vapor | H2O | 80.44 | 69.1 |
| vapor | DME | 2.26 | 2.3 |
| vapor | CO2 | 1.27 | 1.2 |

DME-rich liquid 中 DME/CO2，以及 vapor 中 DME/CO2 的 AARD 几乎逐项复现
公开结果。water-rich liquid 中 H2O 的误差也只有约 2%。高 AARD 集中在
DME-rich liquid 和 vapor 的痕量 H2O，以及 water-rich liquid 的少量 DME/CO2；
这是相对误差分母很小和公开 CPA 参数拟合能力共同造成的，不影响三相逸度相等与
总体物料守恒。

## 判断边界

本算例证明的是：生产 CPA 参数接口能够表达公开的交叉缔合模型，三相 Flash
能够在公开 VLLE 条件下找到满足平衡方程的三个相，并得到与原 CPA 工作相近的
组成精度。它不是“所有 CPA 体系都已验证”的证明，也不验证储层流动离散、井模型
或时间推进。若用于流动算例，还应另做密度、黏度和动态质量守恒检查。

## 参考文献

[1] LAURSEN T, RASMUSSEN P, ANDERSEN S I. VLE and VLLE measurements of
dimethyl ether containing systems[J]. Journal of Chemical & Engineering Data,
2002, 47(2): 198-202. DOI: [10.1021/je010154+](https://doi.org/10.1021/je010154%2B).
[包含原始数据表的 DTU 公开学位论文](https://backend.orbit.dtu.dk/ws/portalfiles/portal/5482520/KT2002-Torben%2BLaursen-Measurments%2Band%2Modeling%2Bof%2BVLLE%2Bat%2BElevated%2BPressures.pdf).

[2] FOLAS G K. Modeling of complex mixtures containing hydrogen bonding
molecules[D]. Lyngby: Technical University of Denmark, 2007. ISBN
978-87-91435-45-5. [DTU 公开全文](https://backend.orbit.dtu.dk/ws/files/123926886/CDocuments_and_SettingsalbDesktopPh.D._Georgios_K._Folas.pdf.pdf).
