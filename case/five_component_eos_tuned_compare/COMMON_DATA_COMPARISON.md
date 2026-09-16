# 三种 EOS 同源公开数据拟合与流动对比报告

## 结论先行

本次比较实现了预定的公平原则：PR、Søreide–Whitson（SW）和 CPA 使用
**完全相同的公开实验记录、训练/验证划分和流动模型**，但分别优化各自
允许调整的参数。

- 含水体系中，SW 的盲验证残差最低：H2O–CO2 为 2.59%，H2O–CH4 为
  11.17%。这正是其水相专用温度/盐度相关式的优势。
- 非水 VLE 与气体密度中，PR/SW 优于当前 CPA 参数包：CO2–C2H6 密度
  误差为 0.606%，三个干气 VLE 体系的合并残差为 6.36%。
- CPA 比普通 PR 更能描述含水体系，但在本数据集上没有超过 SW；它也不
  是所有性质的统一最优模型。
- 在相同 20×20×5 流动算例中，PR 与 SW 的响应几乎重合。CPA 给出明显
  更高的油相饱和度、更低的气相饱和度和更高的生产质量流量；其运行时间
  也比 PR 高约 39%。

因此，不能给出脱离流体类别的单一“冠军”。当前证据支持的结论是：
**SW 最适合本报告覆盖的含盐水–气体系；PR/SW 对干气 VLE 和气体密度更好；
CPA 的流动结果差异主要来自相态分配和水相物性，而不是求解失败。**

## 1. 公平比较方案

### 1.1 固定项与可调项

三种 EOS 共用 201 条公开记录。数据在每条等温线内按固定种子
`20260825` 做确定性 3:1 划分，共 149 条训练记录和 52 条盲验证记录。
验证集不参与参数优化。

| 项目 | PR | SW | CPA |
|---|---:|---:|---:|
| 干气组分对 | 常数 `kij` | 常数 `kij` | 常数 `kij` |
| H2O–CO2、H2O–CH4 | 常数 `kij` | 公开水相相关式的常数偏移 | 常数 `kij` |
| 纯组分参数 | 固定 | 固定 | 固定 |
| 缔合参数 | 不适用 | 不适用 | 固定为公开参数包 |
| 优化器、边界和停止条件 | 相同 | 相同 | 相同 |

这不是让三种模型使用相同的 `kij`，而是让它们在相同证据和相同规则下
各自达到局部最优。未被公开数据覆盖的 H2O–C2H6、H2O–nC4H10、
CO2–nC4H10 和 C2H6–nC4H10 仍使用原参数包先验值，未参与拟合。

### 1.2 实验数据与条件

| 组分体系 | 实验量 | 总数（训练/验证） | 温度 [K] | 压力 [MPa] | 来源 |
|---|---|---:|---:|---:|---|
| CH4–C2H6 | `P-x-y` VLE | 17（12/5） | 203.22–243.61 | 2.124–6.885 | May 等 |
| CH4–nC4H10 | `P-x-y` VLE | 20（13/7） | 203.25–273.42 | 1.311–10.132 | May 等 |
| CO2–CH4 | `P-x-y` VLE | 37（28/9） | 293.128–303.145 | 5.727–7.931 | Petropoulou 等 |
| H2O–CH4 | `P-x-y` VLE | 22（16/6） | 283.89–323.56 | 4.78–19.49 | Frost 等 |
| CO2–C2H6 | 气相 `P-ρ-T-x` | 93（71/22） | 273.15–323.15 | 0.518–6.014 | Yang 等 |
| H2O–CO2–NaCl | CO2 溶解度 | 12（9/3） | 323.15–423.15 | 5.01–20.01 | Messabeb 等 |

所有原始文件来自 NIST ThermoML，下载文件的 DOI 与 SHA-256 均固化在
`reference_data/manifest.json`。盐水数据只取与流动算例一致的
1 mol/kg NaCl 子集。

### 1.3 泡点、露点、相组成和密度如何同时进入拟合

二元 VLE 的每条记录同时包含实验温度 `T`、压力 `P`、液相组成 `x` 和
气相组成 `y`。在固定 `T,P,x,y` 下，对两个组分计算

```text
r_i = ln(f_i^L / f_i^V)
```

理想拟合时两个残差都为零。这里 `P-x` 是泡点支，`P-y` 是露点支；两个
组分残差又同时检验液相和气相组成。因此 VLE 拟合不是只看“某条曲线”，
而是对同一条 tie-line 的泡点、露点与组成做联合约束。密度记录使用
`ln(ρcalc/ρexp)`。图 1 的百分比为 `|exp(r)-1|×100%` 的均值；密度点
则退化为通常的绝对相对误差。该指标应称为**归一化残差**，不是直接求解
泡点压力后得到的压力 AARD。

CO2–NaCl 水溶液论文未报告气相水含量，所以该子集只用 CO2 的等逸度
残差，不能声称验证了水蒸气组成。

## 2. 拟合参数

| 组分对 | PR | SW | CPA |
|---|---:|---:|---:|
| H2O–CO2 | -0.03918 | -0.000238（相关式偏移） | 0.15841 |
| H2O–CH4 | -0.20000 | -0.03784（相关式偏移） | 0.02041 |
| CO2–CH4 | 0.11001 | 0.11001 | 0.00962 |
| CO2–C2H6 | 0.17433 | 0.17433 | 0.06771 |
| CH4–C2H6 | 0.00228 | 0.00228 | -0.00460 |
| CH4–nC4H10 | 0.01088 | 0.01088 | -0.02189 |

PR 的 H2O–CH4 参数触及下限 -0.20，说明普通 PR 加单个常数 `kij` 仍无法
吸收该体系的强非理想性。这个边界值不应移植为通用物性参数。PR 与 SW
的干气参数完全相同，是因为当前 SW 后端在无水相时使用同一个 PR 干气
内核；这也是它们干气结果一致的原因。

## 3. 盲验证结果

| 体系 | PR [%] | SW [%] | CPA [%] | 本组最低 |
|---|---:|---:|---:|---|
| H2O–CO2 | 40.45 | **2.59** | 11.04 | SW |
| H2O–CH4 | 185.03 | **11.17** | 22.74 | SW |
| CO2–CH4 | **10.78** | **10.78** | 14.18 | PR/SW |
| CO2–C2H6 密度 | **0.606** | **0.606** | 0.877 | PR/SW |
| CH4–C2H6 | **1.08** | **1.08** | 1.50 | PR/SW |
| CH4–nC4H10 | **4.43** | **4.43** | 6.63 | PR/SW |

按预测残差数加权后：PR、SW、CPA 的全验证集平均值分别为 33.19%、
5.34% 和 8.71%。这个总数会被 PR 的 H2O–CH4 大误差强烈支配，因此只能
作为数据集级摘要，不能代替分体系判断。CO2–CH4 最大误差集中在临界区，
单个温度无关 `kij` 难以同时描述三条等温线的临界邻域。

图 1 直接展示各体系验证残差，图 2 检查 CO2–C2H6 密度的计算–实验一致性，
图 3 给出干气 `kij`。图均无大标题、无网格、单图窗、四边框，图例位于框内；
中文版本的图例仍为英文。

## 4. 加入流动后的结果

### 4.1 共同计算条件

| 项目 | 设置 |
|---|---|
| 网格/域 | 20×20×5；1000×600×50 m |
| 初始状态 | 60 bar，305 K，`z=[0.30,0.10,0.15,0.15,0.30]` |
| 盐度 | 1 mol/kg NaCl |
| 注入井 | CO2 注入；参考体积流量 100000 m³/day |
| 生产井 | 52 bar BHP |
| 边界 | 外边界无流 |
| 时长 | 5 day；20 个 0.25 day 输出区间 |
| 并行 | MPI 2 ranks |

### 4.2 物理响应与计算代价

| 指标 | PR | SW | CPA |
|---|---:|---:|---:|
| 初始平均油/气/水饱和度 | 0.7175 / 0.1996 / 0.0829 | 0.7156 / 0.2012 / 0.0832 | 0.8239 / 0.0989 / 0.0772 |
| 第 5 天平均油/气/水饱和度 | 0.7160 / 0.2012 / 0.0829 | 0.7140 / 0.2027 / 0.0832 | 0.8219 / 0.1009 / 0.0772 |
| 初始生产质量流量 [kg/s] | 17.280 | 17.202 | 21.913 |
| 第 5 天生产质量流量 [kg/s] | 7.383 | 7.353 | 8.973 |
| SNES / KSP 迭代 | 60 / 328 | 60 / 328 | 60 / 326 |
| 模拟墙钟时间 [s] | 172.23 | 175.06 | 239.15 |
| 接受/拒绝时间步 | 20 / 0 | 20 / 0 | 20 / 0 |
| 最大组分质量守恒相对误差 | 5.04e-14 | 4.19e-14 | 8.82e-14 |

CPA 比 PR 慢 38.85%，但非线性迭代数相同、线性迭代数略少，所以额外耗时
主要来自每次热力学状态评估中的缔合求解，而非流动方程更难收敛。三组
质量守恒均达到约 1e-13，且没有拒绝时间步，说明差异是稳定的模型输出，
不是数值发散伪影。

CPA 的气相饱和度显著较低，说明在当前 305 K、60 bar、五组分总体组成下，
它把更多物质分配到液相。相态比例、相密度和黏度共同改变井的质量流量，
因此不能只用体积产量比较 EOS。

## 5. 适用边界

本实验已经能回答“在当前公开数据覆盖范围内，各 EOS 独立优化后谁拟合
得更好”。它还不能证明某个 EOS 对完整五组分储层流体普遍最优，原因是：

1. 没有找到同时报告本五组分全部二元对、泡/露点、两相组成和密度的单一
   公开数据集；本报告使用可追溯的组合数据集。
2. 流动初始流体含 30 mol% nC4H10，但 CO2–nC4H10、H2O–nC4H10 和
   C2H6–nC4H10 没有进入本轮拟合，CPA 与 PR/SW 的流动差异包含这些先验
   参数的外推影响。
3. 当前拟合只调整每个组分对的一个常数，不包含温度相关 `kij`、体积平移
   或更复杂混合规则。
4. 墙钟时间是同一台机器的一次完整运行，用于量级判断；严谨性能结论需要
   交错顺序的多次重复运行并报告中位数与离散度。

下一轮若要把“流动差异”更严格地归因于 EOS，最有价值的新数据是
305 K、约 50–65 bar 附近的 CO2–nC4H10 与含水五组分相平衡/密度数据。

## 6. 图与可复现文件

- 完整中英文图件为可再生成产物，不再纳入 Git。
- 每张图的绘图数据仍保留于：`figures/common_fit/source_data/`。
- 需要图件时应由对应绘图流程从保留数据重建。
- 冻结实验表：`tools/example/five_component_common_fit/data/common_calibration_data.csv`
- 数据划分与转换：`tools/example/five_component_common_fit/data/dataset_manifest.json`
- 参数与计算入口：`tools/example/five_component_common_fit/main.cpp`

## 参考文献

[1] MAY E F, GUO J Y, OAKLEY J H, et al. Reference quality vapor-liquid equilibrium data for the binary systems methane + ethane, + propane, + butane and, + 2-methylpropane, at temperatures from (203 to 273) K and pressures to 9 MPa[J]. *Journal of Chemical & Engineering Data*, 2015, 60(12): 3606-3620. DOI: 10.1021/acs.jced.5b00610.

[2] PETROPOULOU E, VOUTSAS E, WESTMAN S F, et al. Vapor-liquid equilibrium of the carbon dioxide/methane mixture at three isotherms[J]. *Fluid Phase Equilibria*, 2018, 462: 44-58. DOI: 10.1016/j.fluid.2018.01.011.

[3] YANG X, BEN SOUISSI M A, KLEINRAHM R, et al. Vapour-phase (p, ρ, T, x) behaviour and virial coefficients for the (ethane + carbon dioxide) system[J]. *The Journal of Chemical Thermodynamics*, 2018, 122: 204-213. DOI: 10.1016/j.jct.2018.02.021.

[4] FROST M, KARAKATSANI E, VON SOLMS N, et al. Vapor liquid equilibrium of methane with water and methanol: Measurements and modeling[J]. *Journal of Chemical & Engineering Data*, 2014, 59(4): 961-967. DOI: 10.1021/je400684k.

[5] MESSABEB H, CONTAMINE F, CEZAC P, et al. Experimental measurement of CO2 solubility in aqueous NaCl solution at temperature from 323.15 to 423.15 K and pressure of up to 20 MPa[J]. *Journal of Chemical & Engineering Data*, 2016, 61(10): 3573-3584. DOI: 10.1021/acs.jced.6b00505.

[6] FRENKEL M, CHIRICO R D, DIKY V V, et al. XML-based IUPAC standard for experimental, predicted, and critically evaluated thermodynamic property data storage and capture (ThermoML)[J]. *Pure and Applied Chemistry*, 2006, 78(3): 541-612. DOI: 10.1351/pac200678030541.
