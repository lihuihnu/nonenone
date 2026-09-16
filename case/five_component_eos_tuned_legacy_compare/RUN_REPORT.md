# 三种全组分 EOS 与 Traditional：30 天流动对比报告

## 1. 摘要

本次把统一计算时长由 5 天延长到 30 天。比较对象为 Our PR、Our SW、Our CPA，
以及水不参与相态转化的 Traditional。四组使用相同地质、井控、初始温压、物理
组分总量和 `2 MPI ranks`。

- 四组均到达 `30 day`，完整输出 `120/120` 个 0.25 天目标时刻；
- Our PR 和 Our SW 在约 25 天的相态变化区触发自适应缩步，分别接受/拒绝
  `130/3` 和 `139/6` 个内部步，之后均恢复并完成计算；
- Our CPA 和 Traditional 全程以 0.25 天步长完成，没有拒绝步；
- 四组最终平均压力范围为 `59.85155–59.86404 bar`，总体压力响应接近；
- Our PR、Our SW 与 Traditional 的最终饱和度和井响应接近，Our CPA 保持更高
  油相饱和度和更低气相饱和度；
- 最大逐组分相对质量守恒误差为 `9.05e-14–8.52e-10`。

## 2. 计算条件

| 项目 | 设置 |
|---|---|
| 网格 | `20 × 20 × 5`，2000 个活动单元 |
| 尺寸 | `1000 × 600 × 50 m` |
| 孔隙度 | `0.10–0.24`，体积平均 `0.17896` |
| 渗透率 | Kx `2–250 mD`；Ky `1–160 mD`；Kz `0.05–18 mD` |
| 初始温压 | `305 K`，`60 bar` |
| 组分 | H2O / CO2 / CH4 / C2H6 / nC4H10 |
| 总组成 | `0.30 / 0.10 / 0.15 / 0.15 / 0.30` |
| SW 盐度参数 | `1 mol NaCl / kg H2O`，用于 SW alpha/BIP |
| 注入井 | 纯 CO2，`100000 reference m³/day`，完成层 `k=0..1` |
| 生产井 | `52 bar` BHP，完成层 `k=3..4` |
| 边界与离散 | 封闭外边界；TPFA；相势上风；Backward Euler |
| 时间 | 120 个 0.25 天统一目标区间，总计 30 天；内部步可自适应 |
| 并行 | 每组 `2 MPI ranks` |

Traditional 从总组成中移除 H2O，再对四个干组分归一化并进行 PR 油气闪蒸；水使用
独立守恒方程和固定 `998 kg/m³` 密度。Our PR、SW、CPA 则让五个组分全部参与
O/G/W 三相平衡。

## 3. 密度闭合

Our PR 和 Our SW 在不改变逸度与相平衡条件的前提下，对相摩尔体积采用

\[
v^{\mathrm{tr}}=\frac{ZRT}{P}-\sum_i x_i c_i.
\]

H2O 平移量分别为 `3.252997096e-6 m³/mol` 和
`3.267331277e-6 m³/mol`，共同锚定 IAPWS-IF97 在 `305 K, 6 MPa` 下的纯水
密度 `997.6771185 kg/m³`。SW 水相还使用
`c_CO2,aq = 1.469715763e-6 m³/mol`，其依据是 Garcia 的 CO2 水溶液表观摩尔
体积关联式。上述修正只进入密度和相体积，不进入 Z、逸度、alpha 或 BIP。

当前 NaCl 仍是 SW 相平衡参数，并不是独立守恒组分。因此结果代表 IAPWS 水基准
加溶解 CO2 修正，不等同于显式盐质量守恒的完整盐水密度模型。

## 4. 初始状态

| 模型 | So | Sg | Sw | 水密度 (kg/m³) |
|---|---:|---:|---:|---:|
| Our PR | 0.726704 | 0.202151 | 0.071144 | 998.286 |
| Our SW | 0.724764 | 0.203787 | 0.071449 | 999.198 |
| Our CPA | 0.823858 | 0.098896 | 0.077246 | 1001.657 |
| Traditional | 0.722962 | 0.206182 | 0.070856 | 998.000 |

体积平移改变摩尔体积，因此会影响由相摩尔数换算得到的饱和度，但不改变逸度平衡。
PR、SW 与 Traditional 的初态接近；CPA 的油相饱和度高约 10 个百分点。

## 5. 第 30 天结果

### 5.1 储层与井响应

| 模型 | 平均压力 (bar) | So | Sg | Sw | 水密度 (kg/m³) | 生产量 (m³/s) | 注入井 BHP (bar) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Our PR | 59.85215 | 0.720471 | 0.208388 | 0.071141 | 998.296 | 0.525618 | 62.7953 |
| Our SW | 59.85196 | 0.718557 | 0.209995 | 0.071449 | 999.229 | 0.527099 | 62.8068 |
| Our CPA | 59.86404 | 0.815975 | 0.106780 | 0.077245 | 1001.670 | 0.327031 | 62.5912 |
| Traditional | 59.85155 | 0.716779 | 0.212371 | 0.070850 | 998.000 | 0.533179 | 62.8120 |

Our PR 与 Our SW 的生产量相差 `0.28%`，水密度相差 `0.933 kg/m³`；SW 密度
修正没有造成异常流动偏移。Traditional 生产量比 Our PR 高 `1.44%`。Our CPA
生产量比 Our PR 低 `37.8%`，与其较低气相比例及不同相组成、密度和流度一致。

### 5.2 CO2 分配与水相性质

下表为 2000 个网格在第 30 天的未加权算术平均；`Δxw` 为第 30 天减初始值。

| 模型 | 水相 CO2 (×10⁻³) | Δxw (×10⁶) | 水/油分配比 | 气/油平衡比 | 水黏度 (mPa·s) |
|---|---:|---:|---:|---:|---:|
| Our PR | 1.414 | 13.08 | 0.01026 | 1.5651 | 1.1410 |
| Our SW | 3.273 | 31.99 | 0.02384 | 1.5678 | 1.1410 |
| Our CPA | 2.410 | 25.08 | 0.01703 | 1.4317 | 1.1801 |
| Traditional | 0.000 | 0.00 | 0.00000 | 1.5676 | 0.8000 |

SW 的水相 CO2 平均摩尔分数是 PR 的 `2.31` 倍，CPA 位于两者之间。CPA 的气/油
CO2 平衡比低于其他模型，与更少气相相互印证。Traditional 水相 CO2 恒为零是
模型定义，不是数值失败。

### 5.3 收敛、守恒与耗时

| 模型 | 输出区间 | 接受/拒绝内部步 | SNES | KSP | 最小步长 (day) | 耗时 (s) | 最大相对守恒误差 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Our PR | 120 | 130 / 3 | 410 | 2217 | 1.717e-2 | 746.9 | 2.30e-12 |
| Our SW | 120 | 139 / 6 | 447 | 2329 | 1.998e-4 | 817.3 | 2.15e-12 |
| Our CPA | 120 | 120 / 0 | 360 | 1926 | 2.500e-1 | 1179.9 | 9.05e-14 |
| Traditional | 120 | 120 / 0 | 695 | 4040 | 2.500e-1 | 600.1 | 8.52e-10 |

PR/SW 的缩步集中在约 25 天的相态变化区。SW 达到更小的最小步长，但仅有 6 次
拒绝，随后恢复至统一目标时刻并完成 30 天。CPA 每步保持 3 次 Newton，但缔合
位点和密度根计算使单次残差更昂贵。Traditional 需要更多 Newton/KSP 迭代。
四组并发运行时系统负载不同，因此 wall time 只作本次运行记录，不是严格性能基准。

## 6. 图像解读

### 6.1 整体响应

![平均压力](figures/translated_density_30day_20260826/zh/01_mean_pressure.png)

四条压力曲线整体重合，CPA 略高；30 天末极差为 `0.0125 bar`。

![生产井地面产量](figures/translated_density_30day_20260826/zh/02_producer_rate.png)

Our PR、Our SW 与 Traditional 构成接近的一组；Our CPA 持续较低。

![注入井井底压力](figures/translated_density_30day_20260826/zh/03_injector_bhp.png)

BHP 先短暂升高，再随压力传播下降；约 25 天后曲率轻微变化，与相态变化时段一致。

### 6.2 相态、物性与计算成本

![水相平均密度](figures/translated_density_30day_20260826/zh/04_water_density.png)

最终水密度从 `998.000` 到 `1001.670 kg/m³`，最大相对差约 `0.368%`。

![最终平均饱和度](figures/translated_density_30day_20260826/zh/05_final_saturation.png)

PR、SW、Traditional 的饱和度接近；CPA 明显偏向油相。

![运行时间](figures/translated_density_30day_20260826/zh/06_runtime.png)

图中为本次完整模拟 wall time。CPA 热力学闭合成本最高；图值不代表无资源竞争的
重复性能基准。

![最终水相 CO2](figures/translated_density_30day_20260826/zh/07_aqueous_co2_mole_fraction.png)

SW 给出最高水相 CO2，CPA 次之，PR 最低。

![水相 CO2 增量](figures/translated_density_30day_20260826/zh/08_aqueous_co2_change.png)

扣除初始闪蒸基线后，SW、CPA、PR 的净增量依次降低。

![CO2 水油分配比](figures/translated_density_30day_20260826/zh/09_co2_water_oil_partition.png)

SW 的水/油分配比约为 PR 的 `2.32` 倍，差异集中在含水体系分配。

![CO2 气油平衡比](figures/translated_density_30day_20260826/zh/10_co2_gas_oil_partition.png)

PR、SW、Traditional 约为 `1.57`，CPA 为 `1.43`。

![水相黏度](figures/translated_density_30day_20260826/zh/11_water_viscosity.png)

Our PR/SW 约 `1.141 mPa·s`，CPA 为 `1.180 mPa·s`，Traditional 固定为
`0.800 mPa·s`；流动差异不能只由密度解释。

### 6.3 第 30 天垂向平均场

每个场图对同一 `(i,j)` 柱的 5 层作未加权平均，不插值、不平滑。同一物理量的
四模型使用共享色标，且每幅图片只有一个图窗。

![Our PR 压力变化场](figures/translated_density_30day_20260826/fields/zh/01_pressure_change_our_pr.png)

压力呈注入端升高、生产端降低的偶极结构。

![Our CPA 气相饱和度场](figures/translated_density_30day_20260826/fields/zh/02_gas_saturation_our_cpa.png)

生产井附近气相饱和度升高；CPA 的全场基值显著低于其他模型。

![Our SW 水相 CO2 场](figures/translated_density_30day_20260826/fields/zh/03_aqueous_co2_our_sw.png)

SW 具有最高水相 CO2 背景，并在注入井附近形成局部富集。

![Our CPA 气相 CO2 场](figures/translated_density_30day_20260826/fields/zh/04_gas_co2_our_cpa.png)

气相 CO2 的强变化集中在注入井附近。

![Traditional 水饱和度变化场](figures/translated_density_30day_20260826/fields/zh/05_water_saturation_change_traditional.png)

水饱和度变化为 `10⁻⁴` 量级并集中在井附近，30 天内未出现贯穿井间的水前缘。

仓库只保留正文直接引用的中文场图预览；完整中英文场图属于可再生成产物，
可使用 `plot_field_maps.py` 从对应结果目录重建。

## 7. 可复现性与验证

- Our PR/SW/CPA 原始结果：
  `../five_component_eos_tuned_compare/results/translated_density_30day_20260826/{pr,sw,cpa}`；
- Traditional 原始结果：`results/translated_density_30day_20260826/pr`；
- 仓库内保留正文直接引用的中文对比图和场图预览；
- 完整中英文 PNG/PDF 由 `plot_model_comparison.py` 与 `plot_field_maps.py` 按需重建，不再全部纳入 Git；
- 绘图数据与来源记录：`source_data.csv`、`composition_statistics.csv`、
  `figure_manifest.json` 和 `fields/figure_manifest.json`；
- 绘图不做滤波、平滑、插值或归一化；
- 验证包括终态时间、121 个储层/井记录、内部步计数、2000 单元初末快照、
  非有限值和逐组分质量守恒检查。

本实验比较的是四个完整模型包在一个固定流动场景中的响应，不构成任意温压和
流体体系下的 EOS 普遍优劣排序。

## 8. 参数与数据来源

1. Peng, D.-Y.; Robinson, D. B. A New Two-Constant Equation of State.
   *Industrial & Engineering Chemistry Fundamentals* **1976**, *15*(1), 59–64.
   https://doi.org/10.1021/i160057a011.
2. Søreide, I.; Whitson, C. H. Peng–Robinson Predictions for Hydrocarbons,
   CO2, N2, and H2S with Pure Water and NaCl Brine. *Fluid Phase Equilibria*
   **1992**, *77*, 217–240. https://doi.org/10.1016/0378-3812(92)85105-H.
3. International Association for the Properties of Water and Steam. *Revised
   Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic
   Properties of Water and Steam*, R7-97(2012), revised 2018.
   https://iapws.org/documents/release/IF97-Rev.
4. Garcia, J. E. *Density of Aqueous Solutions of CO2*; LBNL-49023, Lawrence
   Berkeley National Laboratory, 2001. https://escholarship.org/uc/item/6dn022hb.
5. Tsivintzelis, I.; Kontogeorgis, G. M. Modelling Phase Equilibria for Acid
   Gas Mixtures Using the Cubic-Plus-Association Equation of State. *The Journal
   of Supercritical Fluids* **2015**, *104*, 29–39.
   https://doi.org/10.1016/j.supflu.2015.05.015.

公共数据拟合过程与冻结参数见
`../five_component_eos_tuned_compare/COMMON_DATA_COMPARISON.md`。
