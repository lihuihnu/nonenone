# H2O–CO2–nC10 二维 benchmark 详细参数文档

## 1. 文档范围与参数优先级

本文档记录当前统一仓库中 H2O–CO2–nC10 规则二维 benchmark 的正式参数基线，覆盖 Traditional、New-PR、New-SW 和 New-CPA 四种配置。文档依据提交 `483e808` 对应的生产代码、当前运行模板以及正式 0.1 PVI 结果日志编写。

发生不一致时，参数解释按以下优先级处理：

1. 算例生产代码 `case_config.hpp`、`benchmark_common.hpp` 和 `well_config.hpp`；
2. 主程序中的初始化与运行时配置；
3. 正式运行的 `run.log`、`simulation_summary.csv` 和守恒输出；
4. 本文档与数值实验报告；
5. 历史讨论、临时脚本和旧交接材料。

需要特别区分：代码默认运行终点为 0.5 PVI，而当前论文图件使用的是运行时覆盖得到的 0.1 PVI 数据。两者均在本文档中单列，不应混用。

## 2. 四种模型的定义

| 图中名称 | 可执行算例 | 选择器 | 组分数 | 相数 | 每单元主未知量 | 相态描述 |
|---|---|---|---:|---:|---:|---|
| Traditional | `h2o_co2_nc10_2d_traditional` | 无 | 2 | 3 | 7 | CO2–nC10 组分油气体系，加独立水相 |
| New-PR | `h2o_co2_nc10_2d_benchmark` | `-eos pr` | 3 | 3 | 11 | H2O–CO2–nC10 全组分 O/G/W，普通 PR |
| New-SW | `h2o_co2_nc10_2d_benchmark` | `-eos sw` | 3 | 3 | 11 | H2O–CO2–nC10 全组分 O/G/W，Søreide–Whitson |
| New-CPA | `h2o_co2_nc10_2d_benchmark` | `-eos cpa` | 3 | 3 | 11 | H2O–CO2–nC10 全组分 O/G/W，显式 PR-CPA |

Traditional 中水不进入 CO2–nC10 的组分相平衡，因此水相 CO2、油相 H2O 和气相 H2O 按模型定义不表示。三个 New 模型将 H2O 作为普通守恒组分，在统一三相稳定性分析和 flash 中求解相间分配。

四种配置均关闭独立的旧式 aqueous-CO2 dissolution 开关。这不意味着 New 模型禁止 CO2 溶于水；New 模型的水相 CO2 已由全组分 EOS/flash 直接求得。该开关仅对应 Traditional 路径中的额外 Henry 型溶解模块。

## 3. 单位与索引约定

### 3.1 基本单位

| 量 | 内部单位 |
|---|---|
| 长度 | m |
| 时间 | s；配置输入时间步使用 day |
| 压力 | Pa |
| 温度 | K |
| 绝对渗透率 | m²；文档同时给出 mD |
| 体积流量 | m³/s |
| 质量密度 | kg/m³ |
| 摩尔密度 | mol/m³ |
| 动力黏度 | Pa·s |
| 摩尔质量 | kg/mol |

代码中的换算常数为：

$$
1\ \mathrm{bar}=10^5\ \mathrm{Pa},
$$

$$
1\ \mathrm{mD}=9.869233\times10^{-16}\ \mathrm{m^2},
$$

$$
1\ \mathrm{day}=86400\ \mathrm{s},
$$

$$
1\ \mathrm{year}=365.25\ \mathrm{day}.
$$

统一气体常数为 8.31446261815324 J/(mol·K)。

### 3.2 网格索引

算例源代码使用从 0 开始的 i、j、k 索引；论文和人工检查通常使用从 1 开始的索引。`input_index` 是从 0 开始的线性输入顺序。流量符号约定为注入为正、采出为负。

## 4. 几何、网格与边界

| 参数 | 数值 |
|---|---:|
| 计算域尺寸 | 300 m × 100 m × 5 m |
| 逻辑网格 | 60 × 20 × 1 |
| 活跃单元数 | 1200 |
| 单元尺寸 | 5 m × 5 m × 5 m |
| 单元体积 | 125 m³ |
| 总体积 | 150000 m³ |
| 网格类型 | 内存生成的规则 StructuredGrid |
| 外边界 | 无外部传递面，即 no-flow |
| 周期边界 | 无 |

网格尺寸由下式得到：

$$
\Delta x=\frac{300}{60}=5\ \mathrm{m},\qquad
\Delta y=\frac{100}{20}=5\ \mathrm{m},\qquad
\Delta z=\frac{5}{1}=5\ \mathrm{m}.
$$

本算例只有一个水平层，所有相邻单元中心的垂向高程差为零，因此面通量中的重力势差在该几何上为零。当前算例未配置分子扩散、机械弥散或数值弥散模型；组分输运来自 Darcy 对流、相态变化和井源汇。

## 5. 岩石参数与孔隙体积

| 参数 | 数值 | 空间分布 |
|---|---:|---|
| 孔隙度 | 0.20 | 均质 |
| Kx | 100 mD | 均质 |
| Ky | 100 mD | 均质 |
| Kz | 100 mD | 均质 |
| 各向异性 | 1:1:1 | 各向同性 |
| 岩石压缩性 | 未配置 | 孔隙体积乘子保持默认 |

总孔隙体积为：

$$
V_p=L_xL_yL_z\phi
=300\times100\times5\times0.20
=30000\ \mathrm{m^3}.
$$

每个单元的孔隙体积为：

$$
V_{p,\mathrm{cell}}=125\times0.20=25\ \mathrm{m^3}.
$$

运行日志显示的渗透率为 100.000003 mD，这是 mD 与 SI 单位往返换算后的显示舍入，不是空间非均质。

## 6. 初始状态

### 6.1 共同宏观初态

| 参数 | 数值 |
|---|---:|
| 初始压力 | 5.16 MPa = 51.6 bar |
| 温度 | 333.15 K = 60 °C |
| 初始油饱和度 | 0.80 |
| 初始气饱和度 | 0 |
| 初始水饱和度 | 0.20 |
| 初始 CO2 | 无宏观 CO2；仅允许数值痕量 |

压力和温度在全域均匀。模拟为等温过程。

### 6.2 Traditional 初始化

Traditional 直接设置：

$$
S_o=0.80,\qquad S_g=0,\qquad S_w=0.20.
$$

油相与用于数值初始化的气相组成均按 CO2、nC10 顺序取：

$$
\boldsymbol{x}^{o}=\boldsymbol{y}^{g}=(0,1).
$$

水作为独立相初始化，不参与油气相的 H2O–nC10 平衡。

### 6.3 New 模型的平衡初始化

New 模型不是把配置文件中的 `{0.20, 0, 0.80}` 直接写成最终总体摩尔分数。主程序在 5.16 MPa、333.15 K 下对 H2O–nC10 进行两相平衡，并在 H2O 总体摩尔分数 0.50–0.85 的区间内二分搜索 70 次，使 EOS 摩尔密度换算后的水饱和度严格等于 0.20。

最终写入正式初态的总体摩尔组成如下，顺序均为 H2O、CO2、nC10：

| 模型 | H2O | CO2 | nC10 | 得到的 O/G/W 饱和度 |
|---|---:|---:|---:|---|
| New-PR | 0.747021130382 | 0 | 0.252978869618 | 0.8 / 0 / 0.2 |
| New-SW | 0.747083912838 | 0 | 0.252916087162 | 0.8 / 0 / 0.2 |
| New-CPA | 0.747063204494 | 0 | 0.252936795506 | 0.8 / 0 / 0.2 |

这些微小差别来自三种热力学模型对 H2O–nC10 相组成和相摩尔密度的不同预测。此做法保证四模型比较的是同一初始体积饱和度，而不是同一总体摩尔分数造成的不同初始含水体积。

## 7. 组分纯物性参数

New 模型的组分顺序固定为 H2O、CO2、nC10；Traditional 删除 H2O 后按 CO2、nC10 排列。

| 组分 | 临界温度 Tc (K) | 临界压力 Pc (MPa) | 临界摩尔体积 Vc (m³/mol) | 偏心因子 | 摩尔质量 (kg/mol) |
|---|---:|---:|---:|---:|---:|
| H2O | 647.096 | 22.064 | 5.5948 × 10⁻⁵ | 0.3443 | 0.018015268 |
| CO2 | 304.1282 | 7.3773 | 9.4118 × 10⁻⁵ | 0.22394 | 0.0440098 |
| nC10 | 617.70 | 2.103 | 6.10 × 10⁻⁴ | 0.4920 | 0.14228168 |

普通 PR 与 PR-CPA 物理项使用：

$$
\Omega_A=0.45724,\qquad \Omega_B=0.07780,
$$

$$
u=1+\sqrt{2}=2.4142135623730951,
\qquad
w=1-\sqrt{2}=-0.4142135623730951.
$$

`eosModelFlag=5` 对应项目中的 PR78 alpha 分支。New-SW 使用其专用 SW alpha 与 aqueous 参数路径，选择器标志为 1。

## 8. 二元交互系数

### 8.1 共同约束

四种模型的 CO2–nC10 立方 EOS 二元交互系数完全相同：

$$
k_{\mathrm{CO_2,nC10}}=k_{\mathrm{nC10,CO_2}}=0.1141.
$$

这是四模型公平比较的锁定参数，不得在某一模型中单独修改。

### 8.2 Traditional

组分顺序为 CO2、nC10：

$$
\mathbf{k}^{\mathrm{Traditional}}=
\begin{bmatrix}
0 & 0.1141\\
0.1141 & 0
\end{bmatrix}.
$$

### 8.3 New-PR

组分顺序为 H2O、CO2、nC10：

$$
\mathbf{k}^{\mathrm{PR}}=
\begin{bmatrix}
0 & 0.1896 & 0.5000\\
0.1896 & 0 & 0.1141\\
0.5000 & 0.1141 & 0
\end{bmatrix}.
$$

### 8.4 New-SW

New-SW 的非水相矩阵与 New-PR 相同：

$$
\mathbf{k}^{\mathrm{SW,non-aq}}=
\begin{bmatrix}
0 & 0.1896 & 0.5000\\
0.1896 & 0 & 0.1141\\
0.5000 & 0.1141 & 0
\end{bmatrix}.
$$

盐度固定为 0 mol/kg H2O。水相角色使用单独的 water–solute BIP：

$$
k_{\mathrm{H_2O,CO_2}}^{\mathrm{aq}}=-0.06609,
\qquad
k_{\mathrm{H_2O,nC10}}^{\mathrm{aq}}=-0.14229.
$$

CO2–nC10 在水相角色中仍保留共同值 0.1141。SW 水 alpha 在本算例温度与零盐度下由下式计算：

$$
\sqrt{\alpha_w}
=1+0.453\left(1-T_{r,w}\right)
+0.0034\left(T_{r,w}^{-3}-1\right),
$$

$$
T_{r,w}=\frac{333.15}{647.096}=0.5148386020,
$$

$$
\sqrt{\alpha_w}=1.2412933840,
\qquad
\alpha_w=1.5408092651.
$$

### 8.5 New-CPA

New-CPA 明确选择 Peng–Robinson 作为 CPA 的 cubic physical term。该选择是算例级 opt-in，不改变项目中其他 CPA 算例的默认 SRK physical term。

组分顺序为 H2O、CO2、nC10：

$$
\mathbf{k}^{\mathrm{CPA}}=
\begin{bmatrix}
0 & 0.10566 & 0\\
0.10566 & 0 & 0.1141\\
0 & 0.1141 & 0
\end{bmatrix}.
$$

## 9. New-CPA 缔合参数

### 9.1 立方物理项参数

| 组分 | a0 (Pa·m⁶/mol²) | b (m³/mol) | c1 |
|---|---:|---:|---:|
| H2O | 0.15782 | 1.4788 × 10⁻⁵ | 0.6736 |
| CO2 | 0.396304062954024 | 2.66669305019507 × 10⁻⁵ | 0.706477452957888 |
| nC10 | 5.73493521447598 | 1.89999348030598 × 10⁻⁴ | 1.07246071633101 |

CO2 与 nC10 的参数由 PR 纯组分式生成：

$$
a_{0,i}=0.45724\frac{R^2T_{c,i}^2}{P_{c,i}},
\qquad
b_i=0.07780\frac{RT_{c,i}}{P_{c,i}}.
$$

CO2 使用：

$$
c_{1,i}=0.37464+1.54226\omega_i-0.26992\omega_i^2.
$$

nC10 的偏心因子大于 0.49，使用 PR78 高偏心因子分支：

$$
c_{1,i}=0.379642+1.48503\omega_i-0.164423\omega_i^2
+0.016666\omega_i^3.
$$

### 9.2 自缔合与交叉缔合

| 组分 | 缔合能 epsilon (J/mol) | 缔合体积 beta | donor sites | acceptor sites |
|---|---:|---:|---:|---:|
| H2O | 16123.0 | 0.069662 | 2 | 2 |
| CO2 | 8061.5 | 0.15182 | 0 | 1 |
| nC10 | 0 | 0 | 0 | 0 |

H2O 采用 4C 自缔合。CO2 没有 donor site，因此不能自缔合；其正的缔合数组值用于 B2 型 H2O–CO2 交叉缔合。显式交叉参数为：

$$
\epsilon_{\mathrm{H_2O-CO_2}}^{\mathrm{cross}}=8061.5\ \mathrm{J/mol},
$$

$$
\beta_{\mathrm{H_2O-CO_2}}^{\mathrm{cross}}=0.15182.
$$

nC10 不参与缔合。径向分布函数选择 `Simplified`。当体系中没有缔合组分参与时，该 PR-CPA 参数化应退化为相同的普通 PR 物理项。

## 10. 密度与黏度闭合

### 10.1 参考/地面参数

相顺序为 Oil、Gas、Water：

| 量 | Oil | Gas | Water |
|---|---:|---:|---:|
| 参考密度 (kg/m³) | 730.0 | 1.8 | 985.4040020947351 |
| 配置中的回退黏度 (Pa·s) | 2.4 × 10⁻⁴ | 1.86 × 10⁻⁵ | 4.69091 × 10⁻⁴ |

参考密度用于地面/储层体积换算与井输出。它们不是把储层油、气密度固定成常数。配置中的油气黏度值为通用回退数组；本算例的 EOS 油气路径实际使用共同的 LBC 风格黏度闭合。

### 10.2 油气相密度和黏度

油气相由当前 EOS/flash 的相组成与压缩因子计算：

$$
\rho_m=\frac{P}{RTZ},
$$

$$
\rho=\rho_m\sum_i x_iM_i.
$$

Traditional、New-PR、New-SW 与 New-CPA 的 Oil/Gas 均进入同一个 LBC 风格输运物性实现。New-CPA 的物理项为 PR，因此在无缔合 CO2–nC10 子系统、相同组成和相同状态下应与 PR 给出一致的 Z、摩尔密度、质量密度和 LBC 黏度。

### 10.3 New 模型的 IAPWS+Garcia 水相体积

三个 New 模型均开启 IAPWS+Garcia 富水相摩尔体积闭合。纯水密度来自 IAPWS-IF97 Region 1；溶解 CO2 的表观摩尔体积使用 Garcia 形式：

$$
\bar V_{\mathrm{CO_2}}(T)
=\left(37.51-9.585\times10^{-2}t
+8.740\times10^{-4}t^2
-5.044\times10^{-7}t^3\right)\times10^{-6}
\ \mathrm{m^3/mol},
$$

其中 t 为摄氏温度。在 60 °C：

$$
\bar V_{\mathrm{CO_2}}=34.7964496\times10^{-6}\ \mathrm{m^3/mol}.
$$

水相摩尔体积为：

$$
V_m^w=x_{\mathrm{H_2O}}^w\frac{M_w}{\rho_{\mathrm{IAPWS}}(P,T)}
+x_{\mathrm{CO_2}}^w\bar V_{\mathrm{CO_2}}(T).
$$

H2O 摩尔质量取 0.018015268 kg/mol。水相中未被该 H2O–CO2 二元闭合支持的第三组分总摩尔分数上限为 1 × 10⁻⁴；超过上限时程序报错，不做无控制外推。该闭合只修改富水相体积/密度，不修改逸度、Z、稳定性分析或 flash 方程。

### 10.4 McBride–Wright 水相黏度

三个 New 模型的 Water 角色统一采用 McBride–Wright H2O–CO2 关联式：

$$
\ln\left(\frac{\mu_w}{\mathrm{mPa\cdot s}}\right)
=a+b\frac{P}{P_0}
+\frac{c+d(P/P_0)}{T/T_0-1}
+e_1\exp\left[-e_2(T/T_0-1)\right]x_{\mathrm{CO_2}}^w.
$$

参数为：

| 参数 | 数值 |
|---|---:|
| a | −3.705013 |
| b | 0.00289258 |
| c | 3.98950 |
| d | −0.00326 |
| e1 | 65.55968 |
| e2 | 2.46811 |
| T0 | 141.5 K |
| P0 | 1.0 MPa |

在初始压力、温度和无溶解 CO2 时：

$$
\mu_w(5.16\ \mathrm{MPa},333.15\ \mathrm{K},x_{\mathrm{CO_2}}^w=0)
=0.4690906609\ \mathrm{mPa\cdot s}.
$$

未被关联式支持的第三组分总摩尔分数同样限制为 1 × 10⁻⁴。

### 10.5 Traditional 水相物性

Traditional 没有水相 CO2。其水密度使用相同 IAPWS-IF97 Region 1 纯水基准，并通过参考密度定义水体积系数：

$$
B_w(P)=\frac{\rho_{w,\mathrm{ref}}}{\rho_{\mathrm{IAPWS}}(P,T)}.
$$

其水黏度使用 McBride–Wright 公式的纯水部分，即固定取水相 CO2 摩尔分数为零。这样 Traditional 与 New 模型共享水物性的纯水基准，但 Traditional 不引入其模型本身不表示的溶解 CO2 修正。

## 11. 相对渗透率

四模型使用完全相同的平方 Corey 型曲线，无滞回。参数为：

| 参数 | 数值 |
|---|---:|
| 束缚水饱和度 Swc | 0.10 |
| 水驱残余油饱和度 Sorw | 0.15 |
| 临界气饱和度 Sgc | 0.02 |
| 气驱残余油饱和度 Sorg | 0.10 |
| 水端点相对渗透率 | 0.30 |
| 气端点相对渗透率 | 0.80 |
| 油端点相对渗透率 | 0.80 |
| Corey 指数 | 水、气、油均为 2 |

水相：

$$
k_{rw}(S_w)=
\begin{cases}
0, & S_w\le 0.10,\\
0.30\left(\dfrac{S_w-0.10}{0.75}\right)^2,
& 0.10<S_w<0.85,\\
0.30, & S_w\ge 0.85.
\end{cases}
$$

气相：

$$
k_{rg}(S_g)=
\begin{cases}
0, & S_g\le 0.02,\\
0.80\left(\dfrac{S_g-0.02}{0.78}\right)^2,
& 0.02<S_g<0.80,\\
0.80, & S_g\ge 0.80.
\end{cases}
$$

油相：

$$
k_{ro}(S_o)=
\begin{cases}
0, & S_o\le 0.10,\\
0.80\left(\dfrac{S_o-0.10}{0.80}\right)^2,
& 0.10<S_o<0.90,\\
0.80, & S_o\ge 0.90.
\end{cases}
$$

当前算例没有配置毛管压力曲线。若后续引入毛管压力，必须作为新的受控物理变量记录，不能在四模型中的某一个单独启用。

## 12. 井位置、完井与井控

### 12.1 井表

| 参数 | CO2_INJ | PROD |
|---|---|---|
| 类型 | 注入井 | 生产井 |
| 0-based 完井索引 | (0, 9, 0) | (59, 9, 0) |
| 1-based 完井索引 | (1, 10, 1) | (60, 10, 1) |
| input_index | 540 | 599 |
| 完井层数 | 1 | 1 |
| 井半径 | 0.10 m | 0.10 m |
| skin | 0 | 0 |
| 主控制 | 总储层体积速率 | 井底压力 |
| 控制目标 | 9.50642634421 × 10⁻⁵ m³/s | 4.66 MPa |
| 注入相/组分 | Gas / pure CO2 | 不适用 |
| BHP 上限 | 7.16 MPa | 未配置 |
| 初始 BHP | 5.16 MPa | 4.66 MPa |

两井均使用垂直 Peaceman 井指数。初始日志报告两井的总井指数均为 1.352468 × 10⁻¹²；该值由 5 m × 5 m × 5 m 单元、100 mD 水平渗透率、0.10 m 井半径和零 skin 自动计算，不是独立调参项。

### 12.2 注入速率与 PVI

注入速率定义为每年注入 0.10 个孔隙体积：

$$
q_{\mathrm{inj}}
=\frac{0.10V_p}{365.25\times86400}
=9.506426344208685\times10^{-5}\ \mathrm{m^3/s}.
$$

等价日注入量为：

$$
q_{\mathrm{inj}}=8.2135523614\ \mathrm{m^3/day}.
$$

累计注入孔隙体积倍数定义为：

$$
\mathrm{PVI}(t)=\frac{q_{\mathrm{inj}}t}{V_p}.
$$

在注入井维持速率控制、未触发 7.16 MPa BHP 上限时：

$$
t=365.25\ \mathrm{day}\Longleftrightarrow0.1\ \mathrm{PVI},
$$

$$
t=1826.25\ \mathrm{day}=5\ \mathrm{year}
\Longleftrightarrow0.5\ \mathrm{PVI}.
$$

注入井以 `ReservoirTotalRate` 控制，目标量是原位各相体积流量之和，不是地面体积流量。当所需 BHP 达到 7.16 MPa 时允许切换到 BHP 上限控制。生产井始终以 4.66 MPa BHP 控制，未设置额外产量限制。

## 13. 关闭的附加物理

| 物理项 | 状态 | 说明 |
|---|---|---|
| 旧式水相 CO2 溶解模块 | OFF | New 模型由全组分 EOS/flash 表示互溶 |
| 竞争吸附 | OFF | thetaMax 与系数数组均为 0 |
| Land 残余气捕集 | OFF | Land 常数为 0 |
| 相对渗透率滞回 | OFF | 使用单值 Corey 曲线 |
| 盐度 | 0 mol/kg H2O | 纯水基线 |
| 分子扩散 | 未配置 | 不进入本算例 |
| 机械弥散 | 未配置 | 不进入本算例 |
| 毛管压力 | 未配置 | 不进入本算例 |

未启用吸附时仍保留的标准状态 101325 Pa、288.15 K，以及岩石密度 2650 kg/m³，均为通用配置占位，不参与本次计算结果。

## 14. 非线性与时间步参数

### 14.1 代码内默认值

| 参数 | 默认值 |
|---|---:|
| fugacity residual scale | 1.0 |
| variable bounds | OFF |
| SNES stagnation guard | ON |
| stagnation 最少迭代数 | 8 |
| stagnation 窗口 | 5 |
| 最小相对改善 | 1 × 10⁻⁴ |
| 固定输出区间数 | 50 |
| 名义输出间隔 | 36.525 day |
| 自适应内部时间步 | ON |
| 最小内部时间步 | 1 × 10⁻⁵ day = 0.864 s |
| 失败切步因子 | 0.5 |
| 成功增长因子 | 1.25 |
| 困难步收缩因子 | 0.8 |
| easy SNES 阈值 | 6 iterations |
| difficult SNES 阈值 | 14 iterations |
| 最大时间步重试 | 16 |
| 最大井控重求次数 | 8 |

在注入井持续维持目标速率时，默认总时长对应的名义注入量为：

$$
t_{\mathrm{default}}=50\times36.525=1826.25\ \mathrm{day}
\Longleftrightarrow0.5\ \mathrm{PVI}.
$$

若未来长时运行触发注入井 BHP 上限，实际 PVI 必须由 `well_history.csv` 中的接受步注入量积分，不应再用名义速率乘总时间代替。

主程序还启用残差缩放，并把 rate-well residual floor 设置为本算例的实际注入速率，避免通用 1 m³/s 量级基准掩盖小流量井的控制误差。

### 14.2 当前受管运行模板的 PETSc 默认值

当前 `run.sh.in` 生成的运行脚本使用：

| 层级 | 设置 |
|---|---|
| SNES | `newtonls`，basic line search |
| SNES max_it | 50 |
| SNES atol | 1 × 10⁻⁵⁰ |
| SNES rtol | 1 × 10⁻⁸ |
| SNES stol | 主程序未显式覆盖时设为 1 × 10⁻¹² |
| KSP | GMRES，restart 50，right preconditioning |
| KSP atol | 1 × 10⁻⁵⁰ |
| KSP rtol | 1 × 10⁻⁵ |
| PC | ASM restrict，overlap 1 |
| 子域求解 | preonly + ILU(0) |
| zero pivot | 1 × 10⁻⁵⁰ |

当前 WSL 机器配置使用 PETSc 3.22.2，`PETSC_ARCH=test1`，运行器优先系统 `mpiexec`。机器路径保存在不入库的 `config/hpc.local.mk`，不属于物理参数。

## 15. 输出配置

默认每个固定输出时刻写出或打印：

- 井汇总、分相井数据和井历史；
- 组分质量平衡、组分库存和质量总量；
- Newton/非线性求解历史；
- 时间步接受与拒绝记录；
- 井控切换历史；
- 储层诊断与最终求解摘要；
- 最终 `solution_final.csv`。

默认不写每个输出时刻的完整 solution snapshot，也不写 phase-state snapshot。最终解始终按输入网格顺序写出。

New 模型未指定结果目录时分别使用 `./results/new-pr`、`./results/new-sw` 和 `./results/new-cpa`。Traditional 默认目录为 `./results/traditional`。

## 16. 当前正式 0.1 PVI 运行实况

论文图件当前采用 365.25 day 的最终场。该组正式运行均为单 MPI rank、开启自适应时间步，井控和物理参数与上文一致。

| 模型 | 结果标签 | 固定输出区间 | 接受/拒绝内部步 | 非线性求解数 | 接受步 SNES | 接受步 KSP | 最小/最大 dt (day) | 运行时间 (s) |
|---|---|---:|---:|---:|---:|---:|---|---:|
| Traditional | `traditional_phase_role` | 2 | 60 / 13 | 73 | 396 | 396 | 0.0278664 / 11.5253 | 32.46 |
| New-PR | `pr_0p1` | 2 | 38 / 13 | 51 | 199 | 199 | 0.178345 / 28.8131 | 83.00 |
| New-SW | `sw_0p1_phase_order_fix` | 2 | 47 / 13 | 60 | 278 | 278 | 0.0713379 / 23.6037 | 115.90 |
| New-CPA | `cpa_0p1_optimized4` | 10 | 65 / 14 | 79 | 387 | 13026 | 0.0713379 / 17.1740 | 758.78 |

最终最大组分相对质量平衡误差为：

| 模型 | 最大相对误差 |
|---|---:|
| Traditional | 6.521 × 10⁻¹⁰ |
| New-PR | 4.399 × 10⁻¹⁰ |
| New-SW | 6.888 × 10⁻¹⁴ |
| New-CPA | 5.709 × 10⁻¹¹ |

四组运行的井控切换文件均只有表头，即井控切换次数为 0。

Traditional、New-PR 与 New-SW 的正式日志明确记录为 `preonly/lu` 串行线性求解。`cpa_0p1_optimized4` 保存了完整摘要和结果，但缺少对应 `run.log`；摘要中的 KSP 迭代数证明它使用了迭代线性求解，但不能仅凭摘要可靠还原全部 PETSc flags。因此，最终场和守恒结论可复核，而若需要逐位复现该 CPA 运行的代数迭代历史，应使用当前受管 GMRES/ASM 模板重新运行并保存完整日志。

另一个必须说明的差异是：前三个模型保存 2 个固定输出区间，New-CPA 保存 10 个。四模型的 0.1 PVI 最终场可直接比较，但现有文件不适合做同等时间分辨率的瞬态曲线比较。

## 17. 推荐复现命令

先按目标机器建立不入库的配置并检查：

```bash
cp config/wsl.example.mk config/hpc.local.mk
make print-config
make doctor-local
```

编译两个算例：

```bash
make case CASE=h2o_co2_nc10_2d_traditional -j
make case CASE=h2o_co2_nc10_2d_benchmark -j
```

复现 0.1 PVI 最终状态时，建议四模型统一使用 10 个输出区间，每个 36.525 day：

```bash
make run CASE=h2o_co2_nc10_2d_traditional NP=1 \
  RESULT_DIR=./results/traditional_0p1 \
  RUN_ARGS='-numSteps 10 -dt 36.525 -adaptive_dt true'

make run CASE=h2o_co2_nc10_2d_benchmark NP=1 EOS=pr \
  RESULT_DIR=./results/new_pr_0p1 \
  RUN_ARGS='-numSteps 10 -dt 36.525 -adaptive_dt true'

make run CASE=h2o_co2_nc10_2d_benchmark NP=1 EOS=sw \
  RESULT_DIR=./results/new_sw_0p1 \
  RUN_ARGS='-numSteps 10 -dt 36.525 -adaptive_dt true'

make run CASE=h2o_co2_nc10_2d_benchmark NP=1 EOS=cpa \
  RESULT_DIR=./results/new_cpa_0p1 \
  RUN_ARGS='-numSteps 10 -dt 36.525 -adaptive_dt true'
```

运行到代码默认的 0.5 PVI 时，不需要覆盖步数和输出间隔：

```bash
make run CASE=h2o_co2_nc10_2d_benchmark NP=1 EOS=cpa \
  RESULT_DIR=./results/new_cpa_0p5
```

在本地 Linux 或超算上只应修改 `config/hpc.local.mk` 中的机器环境、启动器、模块和路径，不应为了适配机器而修改本文档中的物理参数。

## 18. 修改参数时的控制变量规则

后续开发或敏感性研究至少应遵守以下规则：

1. 四模型必须共享网格、孔渗、初始压力/温度、目标饱和度、相渗和井控；
2. CO2–nC10 的 0.1141 必须保持共同，除非建立一组明确标注的新敏感性算例；
3. IAPWS+Garcia 密度与 McBride–Wright 黏度是锁定的共同水物性闭合；
4. PR-CPA 只能显式 opt-in，不能静默替换其他 CPA 算例的 SRK 默认；
5. 改变一个 EOS 参数时，必须记录矩阵顺序、单位、适用相角色和温度；
6. 四模型瞬态对比必须使用相同输出时刻；
7. 结果比较前必须核对质量守恒、井控切换和实际累计 PVI；
8. 任何毛管压力、重力层差、扩散、弥散、吸附或捕集扩展都应作为新的实验因素，不应混入当前基线。

## 19. 代码溯源

| 内容 | 权威文件 |
|---|---|
| 网格、岩石、PVI、井控、相渗 | `case/h2o_co2_nc10_2d_benchmark/benchmark_common.hpp` |
| New-PR/SW/CPA 组分与 EOS 参数 | `case/h2o_co2_nc10_2d_benchmark/case_config.hpp` |
| New 模型初始化与选择器 | `case/h2o_co2_nc10_2d_benchmark/h2o_co2_nc10_2d_benchmark.cpp` |
| New 模型井组分配置 | `case/h2o_co2_nc10_2d_benchmark/well_config.hpp` |
| Traditional 参数 | `case/h2o_co2_nc10_2d_traditional/case_config.hpp` |
| Traditional 初始化与水物性 | `case/h2o_co2_nc10_2d_traditional/h2o_co2_nc10_2d_traditional.cpp` |
| IAPWS+Garcia 实现 | `models/include/natural/thermo/aqueous_volume.hpp` |
| McBride–Wright 实现 | `models/include/natural/properties/aqueous_viscosity.hpp` |
| LBC 风格油气黏度 | `models/include/natural/properties/compositional_properties.hpp` |
| 当前 PETSc 运行模板 | `case/run.sh.in` 与 `case/Makefile` |
| 0.1 PVI 结果解释 | `docs/H2O_CO2_NC10_0P1_NUMERICAL_EXPERIMENT.md` |

本文档是参数台账，不替代 `docs/MODEL.md` 中的控制方程说明，也不替代 `docs/THERMODYNAMICS.md` 中的 flash、相态切换和后备求解流程说明。
