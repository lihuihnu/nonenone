# 超临界水–干酪根裂解产物二维算例与物性标定设计

## 1. 当前结论

当前程序已经增加两套共享几何和井控的 `60×20×1` 均质二维算例：

| 算例 | 组分 | 可切换热力学模型 | 当前用途 |
|---|---|---|---|
| `scw_kerogen_squalane_2d` | H2O / squalane | PR / SW / CPA | 单重质拟组分、相平衡与黏度机制基线 |
| `scw_kerogen_lmh_2d` | H2O / nC4 / nC10 / squalane | PR / SW / CPA | 轻–中–重组分分异与产出组成 |

这两套算例已经完成编译、初态 PR/SW/CPA 闪蒸预检和 PR 一步流动预检。它们现在适合开展模型筛选和数值敏感性研究，但还不能称为“完成实验标定的预测模型”。主要原因有三点：

1. 目标温压窗内的 H2O–squalane 数据只有四条非临界共存 tie-line，且经典 PR 的单一对称 `kij` 不能同时把两个共存端点压入原文精度带；
2. H2O–nC4 与 H2O–nC10 已找到可用的 Søreide–Whitson、PR 和 CPA 文献先验，但逐点原始数据尚未整理成可回归表；现有参数不能跨 EOS、跨相角色直接移用；
3. 当前黏度仍为 LBC 型预测，squalane 临界体积等输入含筛选性质；文献检索支持把水相改为 IAPWS 2008，并用 f-theory/Pedersen 对烃相建立主模型、把 LBC 降为对照模型。

因此，本轮最重要的改进不是给参数贴上“已标定”的标签，而是把**实验数据层、热力学回归层、流动预测层和模型结构不确定性层分开**。这样能够得到更物理、也更容易解释的现象。

## 2. 二维几何、岩石与井位

| 项目 | 数值 |
|---|---:|
| 网格 | `60×20×1`，1200 个活动网格 |
| 尺寸 | 1.20×0.10×0.10 m |
| 单元尺寸 | 0.020×0.005×0.10 m |
| 孔隙度 | 0.35，均质 |
| 水平渗透率 | `Kx=Ky=1500 mD` |
| 垂向渗透率 | `Kz=150 mD`；单层算例中不产生垂向运移 |
| 孔隙体积 | 0.0042 m³ |
| 注采速率 | 各 0.25 PV/day，对等地层体积速率 |
| 基准温压 | 653.2 K，27.74 MPa |
| 时间 | 40 个 0.10 day 输出区间，共 1.0 PVI |

注入井位于左侧短边 `i=0,j=9,k=0`，生产井位于右侧短边 `i=59,j=9,k=0`。`ny=20` 为偶数，短边几何中心位于 `j=9` 和 `j=10` 两个网格中心之间；主算例使用 `j=9`，相对中心线偏移 2.5 mm。严格的数值对称性验证应增加以下两项，而不是把当前井位称为“无偏差中心井”：

- 镜像井位 `j=10`，比较采出曲线和压力差；
- `60×21×1` 奇数横向网格，将井严格放在 `j=10`。

二维网格保持了一维算例的外形、孔隙体积和 PVI 速率，因而一维结果可作为“无横向扩展”的控制组。井半径取 1 mm，以满足细化后 5 mm 横向单元中的 Peaceman 井指数几何约束。

## 3. H2O–squalane PR 标定与验收

### 3.1 进入目标窗的数据

Stevenson 等（1994）的目标窗为 637.2–653.2 K、26.5–29.5 MPa。当前已逐表录入四条非临界 LLE tie-line，组成均为水摩尔分数：

| T / K | P / MPa | 水富相 `xH2O` | squalane 富相 `xH2O` | 用途 |
|---:|---:|---:|---:|---|
| 637.2 | 26.88 | 0.9997 | 0.809 | 标定 |
| 653.2 | 27.05 | 0.9938 | 0.932 | 验证 |
| 653.2 | 27.74 | 0.9966 | 0.911 | 标定 |
| 653.2 | 29.47 | 0.9983 | 0.885 | 验证 |

原文给出的组成测量精度为水富相优于 ±0.0005、squalane 富相 ±0.002，但没有给出覆盖因子和置信水平。因此这些数值在当前工作中只作为确定性验收带，不能直接解释成 `1σ` 标准差。642.2 K 的数据表没有落入上述压力窗的非临界 tie-line，不能人为补出第三个温度层。

### 3.2 建议回归式和目标函数

对经典 PR 先只允许两个可解释参数：

```text
kij(T) = kref + b * (1/T - 1/653.2 K)
```

目标函数分别标准化两个共存端点，而不是只拟合其中一个相：

```text
J = sum_cal [((xw,WR,calc-xw,WR,exp)/0.0005)^2
           + ((xw,SR,calc-xw,SR,exp)/0.002)^2]
```

用 637.2 K、26.88 MPa 和 653.2 K、27.74 MPa 两条 tie-line 标定，用 653.2 K 的另外两条压力点验证。必须同时报告端点残差、相数、密度根和失败闪蒸，不得只报告优化器的 `J` 值。

当前离线检查已经表明：单一经典 PR `kij` 无法同时满足两支组成验收带。这是**模型结构失配**，不能靠继续增加小数位解决。当前流动基线使用组成目标直接拟合的 `kij(653.2 K)=0.0532336595` 和逆温斜率 `1005.04545 K`。它改善了 squalane 富相端点，但水富相端点和端点逸度闭合仍未通过实验验收，因此必须与已发表的 `0.2395` 参数化共同构成敏感性带，不能把组成目标的最小二乘收敛写成完整热力学标定成功。

## 4. H2O–nC4 与 H2O–nC10 温度相关 BIP

### 4.1 当前实现

当前多组分算例在 653.2 K 使用以下非水相 PR 锚点：

| 组分对 | `kij(653.2 K)` | 温度形状 | 状态 |
|---|---:|---|---|
| H2O–nC4 | 0.5091 | 常数 | Søreide–Whitson 非水相公开先验，尚未直接回归 |
| H2O–nC10 | 0.2618373654 | `-0.1646700201+278.5946242/T` | 573.2–613.2 K 文献点回归后向目标窗外推 |
| H2O–squalane | 0.0532336595 | `kref+1005.04545(1/T-1/653.2)` | 目标窗组成拟合；模型结构失配仍在 |

一个重要的数值发现是：目标响应对 BIP 不是平滑弱敏感，而是会跨越相边界。当前所有模型统一采用相同总体组成，并要求初态为 Oil 角色单相；水相只能由注入后的稳定性判据生成。因此 BIP 必须先在二元闪蒸中回归，再锁定后进入多组分流动，不能在流动历史拟合中让它任意补偿相对渗透率或黏度误差。

### 4.2 可直接复用的公开 BIP 与适用边界

以下参数已在 2026-09-16 从原论文或出版商页面复核。最重要的限定是：**BIP 属于“EOS + 混合规则 + 相角色 + 纯组分参数”这一整套模型，不能只抄一个数字。**

Søreide–Whitson PR 对水–烃采用两套参数。非水相使用常数：

| 烃组分 | `kij(non-aq)` |
|---|---:|
| CH4 | 0.4855 |
| C2H6 | 0.4920 |
| C3H8 | 0.5525 |
| nC4 | 0.5091 |
| C5+ 未单列烃 | 0.5000 |

因此当前算例的 nC10 `0.5000` 有明确的 Søreide–Whitson **非水相默认值**依据；nC4 应从 `0.5000` 改为原表的 `0.5091`。但水富相不能继续使用这些数值，而应使用 Søreide–Whitson 的温度/盐度关联。零盐度下，按当前纯组分 `Tc` 和 `omega` 计算得到：

| 组分 | `kij,aq(637.2 K)` | `kij,aq(653.2 K)` | 状态 |
|---|---:|---:|---|
| nC4 | 0.123783 | 0.128315 | 原关联组分范围内；目标温度仍需用 Tian 数据验证 |
| nC10 | 0.063495 | 0.065009 | 超过原关联的 nC4 组分上限，只能作为外推先验 |
| nC16 | 0.083184 | 0.082788 | 更重烃外推，不作为标定结果 |

这两套数值相差很大并不矛盾：它们分别服务于水富相和非水相的逸度计算。当前普通 PR backend 对所有相使用同一个对称 `kij`，所以无法忠实复现这套文献模型；只有 SW 路径已经具备水相专用 BIP。若希望普通 PR 也进行公平比较，需要增加按 thermodynamic role 选择 BIP 的能力，或明确把它定义成“单一 BIP 的经典 PR”结构对照。

对本研究最有用的直接 PR 回归值来自 Teratani 等：

| 二元体系 | T / K | `kij` (`lij=0`) | 说明 |
|---|---:|---:|---|
| H2O–nC10 | 573.2 | 0.3157 | 标准 PR、`lij=0` 的温度逐点回归，低于目标窗 |
| H2O–nC10 | 593.2 | 0.3167 | 标准 PR、`lij=0` 的温度逐点回归，低于目标窗 |
| H2O–nC10 | 613.2 | 0.2836 | 标准 PR、`lij=0` 的温度逐点回归，低于目标窗 |
| H2O–squalane | 637.2 | 0.2207 | 与目标窗直接重合 |
| H2O–squalane | 642.2 | 0.2286 | 与目标温区直接重合 |
| H2O–squalane | 653.2 | 0.2395 | 当前算例已采用 |

nC10 三点说明单一经典 PR 的有效 `kij` 具有模型依赖和温度依赖；按 `a+b/T` 回归在 653.2 K 得到 `0.2618373654`，只能作为目标窗外推值，不能当作测量值。squalane 的已发表参数点与当前直接组成拟合给出明显不同的有效参数，正是 PR 结构不确定性带的一部分。

Oliveira 等的 SRK-CPA 给出了另一套可复用参数：

| n-烷烃 | CH4 | C2 | C3 | nC4 | nC5 | nC6 | nC7 | nC8 | nC9 | nC10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `kij` | 0.165 | 0.141 | 0.117 | 0.092 | 0.068 | 0.044 | 0.019 | -0.005 | -0.0293* | -0.054 |

`*` nC9 为论文关联 `kij=-0.0243 Cn+0.1894` 的插值；其余为论文表值。该关联只由 CH4–nC10 数据建立，且物理项是 **SRK**。当前算例选择的是 PR-CPA，因而不能直接粘贴这些数值；可执行方案是增加一个 SRK-CPA 文献基线，或在当前 PR-CPA 中重新回归。将关联外推到 nC16 和 C30 会分别得到约 -0.199 和 -0.540，已经远离原标定域，不应使用。

据此，当前三模型最科学的参数中心为：

| 模型 | H2O–nC4 | H2O–nC10 | H2O–squalane |
|---|---|---|---|
| 经典 PR 非水相 | 0.5091 | `-0.1646700201+278.5946242/T` | 当前组成拟合关系；另保留已发表参数化作结构敏感性 |
| SW 水相 | 原关联 | 关联外推，仅作敏感性 | 不使用 SW 绝对值；仅作结构外推带 |
| SRK-CPA 文献基线 | 0.092 | -0.054 | 无直接参数，必须回归 |
| 当前 PR-CPA | 不移用 SRK-CPA 数值 | 不移用 SRK-CPA 数值 | 用目标 tie-line 独立回归 |

### 4.3 完成直接回归所需的数据动作

H2O–nC4 文献覆盖到 695 K 和 306 MPa，可支持目标温区附近的直接回归；H2O–nC10 的公开 LLE/VLE 数据主要在 573.2–613.2 K、压力至约 23 MPa，向 637.2–653.2 K 外推不可避免。下一步应：

1. 从原始表格录入 `T/P/总体组成/两相组成/相型/测量不确定度`，保留来源表号和点号；
2. 每个二元体系使用 `kij(T)=a+b/T`，如残差对压力有系统趋势，再比较 `a+b/T+cP`，但不在同一轮同时加入更多无物理约束参数；
3. 按温度整组留出，而不是随机拆散相邻压力点：最低或最高温度整组作外推验证；
4. 对 nC10 报告“目标窗外推带”，其宽度至少包含回归协方差、数据不确定度和 PR/SW/CPA 结构差异；
5. 将“预测相数错误”设为硬失败，不能以组成残差数值惩罚替代。

## 5. 高温高压密度与黏度校准

当前生产路径使用 LBC 型相黏度。它能生成随相组成、密度和温压变化的黏度，但 squalane 的临界体积仍含筛选性质，因此当前的绝对油相黏度只能用于机制筛选。

### 5.1 密度：IAPWS + 一致体积平移

文献使用最广、且与当前程序能力最匹配的分层方案是：

1. **水富相**：纯水基准使用 IAPWS-95。其有效域覆盖本算例；不要用调整 PR 的 `Tc/Pc` 去同时补偿水密度和相平衡。
2. **烃富相**：PR/SRK/CPA 的液体密度采用 Péneloux 型组分体积平移，`v_translated = v_EOS - sum(xi ci)`。常数平移不改变相平衡，因此 `ci` 只对纯组分/单相密度回归并在 BIP 标定时锁定。
3. **重质拟组分**：可用 Jhaveri–Youngren 的三参数 PR/volume-shift 思路生成 C7+ 初值，再用本体系密度数据微调；不能让体积平移去吸收错误相数。
4. **squalane 数据边界**：Mylona 等的高精度 Tait 密度关联覆盖 273–473 K、最高 200 MPa，可用于检查常温到中温参数，但 637–653 K 是明确外推，不能据此宣称超临界水温区密度已验证。

当前 cubic EOS 已实现“只修正密度、不改变 fugacity/Z 根”的 Péneloux 风格路径，所需工作是为 nC4、nC10、squalane 分别回归 `ci`，而不是改写相平衡参数。水富相则继续由 IAPWS 密度作为权威基准，避免重复修正。

### 5.2 黏度：IAPWS + f-theory/Pedersen，LBC 作为对照

建议采用三级校准，严格避免把相平衡误差吸收到黏度参数中：

1. **纯水层**：用 IAPWS 2008 水黏度公式覆盖 637.2–653.2 K、26.5–29.5 MPa；该公式覆盖稳定流体区至 1173 K、300 MPa。水富相不再用 LBC 拟合。
2. **纯烃层**：主模型采用与 cubic EOS 压力项耦合的 f-theory；它已在原始论文中用于 CH4–nC10 及其混合物，并扩展到原油。Pedersen 对应态模型作为独立工程对照；当前 LBC 保留为第三条模型带，而不是唯一预测器。
3. **混合物层**：相组成仍由 PR/SW/CPA flash 给出，黏度模型在各自相组成上计算；不得在水黏度和油黏度之间做线性插值。用至少两个组成、三个压力、三个温度的数据检验混合规则，并留出一个完整温度层验证。
4. **squalane 数据边界**：Mylona 等的两套参考黏度关联覆盖最高 473 K、200 MPa；可用于确定低温参数和检查趋势，但对 637–653 K 只能算外推先验。当前零流动标定中 LBC 对 squalane 的训练/验证 AARD 约为 80%/101%，已经足以否定其作为本算例绝对黏度主模型。

黏度目标函数建议在对数空间构造：

```text
Jmu = sum [log(mu_calc/mu_exp) / sigma_log_mu]^2
```

这样不会让高黏度点仅因数值更大而支配回归。验收至少报告 AARD、最大相对误差、按相/温度分层的偏差和外推标记。没有混合物黏度数据时，应将纯组分校准和混合规则不确定性分开报告，不应声称“混合物黏度已校准”。

## 6. PR/SW/CPA 模型差异带

每一个完全相同的 `P/T/z`、网格、井控和数值容差条件分别运行 PR、SW、CPA。对任一响应 `y(t)` 定义结构差异带：

```text
y_low(t)  = min(y_PR, y_SW, y_CPA)
y_high(t) = max(y_PR, y_SW, y_CPA)
y_mid(t)  = median(y_PR, y_SW, y_CPA)
```

重点响应包括：

- 生产井水质量分数达到 10% 和 50% 的 PVI；
- 1 PVI 的总烃采收率；
- 生产井轻/中/重烃质量分数及轻重选择性；
- 油相、水相最小/平均/生产井黏度；
- 注采压差、最大井底压力和井控切换；
- 水前缘宽度、横向扫掠面积和相态单元比例；
- 每一组分最大相对质量守恒误差。

该 min–max 范围是**模型结构包络**，不是统计置信区间。只有在实验误差、参数协方差和模型权重均有依据时，才能进一步构造概率区间。SW 对 squalane、nC10 存在明确外推，CPA 当前只包含标准水缔合项且未重新拟合重烃交叉项，因此模型差异带尤其需要保留这些适用域标签。

## 7. 优化后的分阶段实验设计

### 阶段 A：无流动 PVT/黏度标定

先完成二元闪蒸和黏度数据回归；流动参数不参与。每套模型保存标定集、验证集、残差和失效相数。只有满足预先规定的相数与残差门槛，才进入阶段 B。

### 阶段 B：确定性数值审计

对基准条件至少计算：

- `60×20,j=9` 主网格；
- `60×20,j=10` 镜像井位；
- `60×21,j=10` 严格中心井；
- `120×40` 或 `120×41` 加密网格。

验收标准建议为：PVI50 变化小于 0.02 PVI、1 PVI 采收率变化小于 2 个百分点、主要组分产出分数变化小于 0.02、质量守恒误差低于 `1e-7`。若未通过，先处理离散误差，不解释 EOS 差异。

### 阶段 C：五因素 Resolution-V 筛选

采用 `2^(5-1)` 半分数设计，生成元 `E=ABCD`，共 16 个角点加 1 个中心曲率检查：

| 因素 | 低水平 | 高水平 | 物理含义 |
|---|---:|---:|---|
| A：温度 | 637.2 K | 653.2 K | 接近/越过水临界温度的相态与黏度变化 |
| B：压力 | 26.5 MPa | 29.5 MPa | 目标实验压力窗 |
| C：注入速率 | 0.10 PV/day | 0.50 PV/day | 平衡接近程度与压降 |
| D：初始水饱和度 | 0.05 | 0.30 | 初始相连通性 |
| E：轻烃质量分数 | 0.05 | 0.25 | 轻组分优先采出与黏度降低 |

每个条件运行 PR/SW/CPA，共 `17×3=51` 个主运行。中心点在确定性模拟中不是“重复测量”，只用于曲率检查；若需要估计实验噪声，应由实验重复提供，而不是重复执行完全相同的确定性求解器。

该 Resolution-V 设计中主效应不与二因子交互混杂，适合先筛选温度、压力、速率、初始水相和轻组分比例。执行顺序已用固定种子随机化并保存为 `flow_screening_design.csv`。当前程序已经能切换 EOS，但批量改变温压、速率、目标初始水饱和度和轻组分比例仍需增加运行时参数及 `P/T/S/z` 反算初始化；因此这 51 个运行是可执行设计，不是已完成结果。

### 阶段 D：局部响应面

只保留阶段 C 中效应最大的三个连续因素，使用中心复合或 Box–Behnken 设计拟合二次响应面。此阶段重点解析 PVI50、采收率和轻/重选择性，不再把所有五个因素同时加密。若模型差异带宽度大于因素效应，则优先补 PVT/黏度实验，而不是增加流动模拟数。

## 8. 当前二维 PR 预检结果

以下仅是 653.2 K、27.74 MPa、0.001 day 的一步流动预检，不代表 1 PVI 预测：

| 指标 | H2O–squalane | H2O–nC4–nC10–squalane |
|---|---:|---:|
| 初始 `beta(O/W)` | 0.89696 / 0.10304 | 0.95782 / 0.04218 |
| 初始 `S(O/W)` | 0.98269 / 0.01731 | 0.99131 / 0.00869 |
| 初始平均油黏度 / Pa·s | `9.361e-4` | `3.253e-4` |
| 初始平均水黏度 / Pa·s | `6.617e-5` | `6.445e-5` |
| 接受/拒绝步 | 1 / 0 | 1 / 0 |
| SNES / KSP 迭代 | 7 / 3213 | 7 / 2057 |
| 最大组分相对质量误差 | `1.28e-13` | `8.02e-10` |

预检说明二维装配、井源项、相平衡、黏度和组分守恒已经耦合工作。此时流体前缘只移动了极小距离，尚不能评价突破时间和生产井组分演化；正式比较必须完成共同的 1 PVI 终点并通过质量守恒门槛。

## 9. 推荐执行顺序

```bash
# 编译二维算例
make case CASE=scw_kerogen_squalane_2d -j2
make case CASE=scw_kerogen_lmh_2d -j2

# 在完全相同条件下形成三 EOS 结构包络
make run CASE=scw_kerogen_squalane_2d NP=2 EOS=pr  RESULT_DIR=./results/pr
make run CASE=scw_kerogen_squalane_2d NP=2 EOS=sw  RESULT_DIR=./results/sw
make run CASE=scw_kerogen_squalane_2d NP=2 EOS=cpa RESULT_DIR=./results/cpa

make run CASE=scw_kerogen_lmh_2d NP=2 EOS=pr  RESULT_DIR=./results/pr
make run CASE=scw_kerogen_lmh_2d NP=2 EOS=sw  RESULT_DIR=./results/sw
make run CASE=scw_kerogen_lmh_2d NP=2 EOS=cpa RESULT_DIR=./results/cpa
```

完整生产运行之前，先执行单元测试和 0.001 day 预检；随后再做 1 PVI。若 PR/SW/CPA 中任一模型出现错误相数、闪蒸失败或质量守恒超限，该点应标为模型失败并保留在包络报告中，不能静默删除。

## 10. 数据与实现入口

- 二维共用网格和井：`case/scw_kerogen_common/benchmark_2d_common.hpp`
- squalane 二维模型：`case/scw_kerogen_squalane_2d/`
- 轻/中/重二维模型：`case/scw_kerogen_lmh_2d/`
- squalane 目标窗数据：`case/scw_kerogen_common/calibration/squalane_target_window.csv`
- 五因素筛选矩阵：`case/scw_kerogen_common/calibration/flow_screening_design.csv`
- 二维装配测试：`test/src/unit/scw_kerogen_squalane_2d_case_test.cpp`、`test/src/unit/scw_kerogen_lmh_2d_case_test.cpp`

## 11. 主要文献

1. Stevenson, R. L., Labracio, D. S., Beaton, T. A., Thies, M. C. (1994). *Fluid-phase equilibria and critical phenomena for the dodecane-water and squalane-water systems at elevated temperatures and pressures*. Fluid Phase Equilibria, 93, 317–336. DOI: `10.1016/0378-3812(94)87016-0`.
2. Søreide, I., Whitson, C. H. (1992). *Peng–Robinson predictions for hydrocarbons, CO2, N2 and H2S with pure water and NaCl brine*. Fluid Phase Equilibria, 77, 217–240. DOI: `10.1016/0378-3812(92)85105-H`.
3. Tian, Y., Zhao, X., Chen, L., Zhu, H., Fu, H. (2004). *High pressure phase equilibria and critical phenomena of water + iso-butane and water + n-butane systems to 695 K and 306 MPa*. Journal of Supercritical Fluids, 30, 145–153. DOI: `10.1016/j.supflu.2003.09.002`.
4. Wang, Q., Chao, K. C. (1990). *Vapor–liquid and liquid–liquid equilibria and critical states of water + n-decane mixtures*. Fluid Phase Equilibria, 59, 207–215. DOI: `10.1016/0378-3812(90)85035-9`.
5. IAPWS (2008, editorial revision 2012). *Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance*, R12-08.
6. Teratani, S., Ota, M., Sato, Y., Inomata, H. (2017). *Development of Predictive Methods of Water–Heavy Oil Phase Equilibrium for Supercritical Water Upgrading Process*. Journal of the Japan Petroleum Institute, 60, 26–33. DOI: `10.1627/jpi.60.26`.
7. Oliveira, M. B., Coutinho, J. A. P., Queimada, A. J. (2007). *Mutual solubilities of hydrocarbons and water with the CPA EoS*. Fluid Phase Equilibria, 258, 58–66. DOI: `10.1016/j.fluid.2007.05.023`.
8. Péneloux, A., Rauzy, E., Fréze, R. (1982). *A consistent correction for Redlich–Kwong–Soave volumes*. Fluid Phase Equilibria, 8, 7–23. DOI: `10.1016/0378-3812(82)80002-2`.
9. Jhaveri, B. S., Youngren, G. K. (1988). *Three-Parameter Modification of the Peng–Robinson Equation of State To Improve Volumetric Predictions*. SPE Reservoir Engineering, 3, 1033–1040. DOI: `10.2118/13118-PA`.
10. Lohrenz, J., Bray, B. G., Clark, C. R. (1964). *Calculating Viscosities of Reservoir Fluids From Their Compositions*. Journal of Petroleum Technology, 16, 1171–1176. DOI: `10.2118/915-PA`.
11. Pedersen, K. S., Fredenslund, A. (1987). *An improved corresponding states model for the prediction of oil and gas viscosities and thermal conductivities*. Chemical Engineering Science, 42, 182–186. DOI: `10.1016/0009-2509(87)80225-7`.
12. Quiñones-Cisneros, S. E., Zéberg-Mikkelsen, C. K., Stenby, E. H. (2000). *The friction theory (f-theory) for viscosity modeling*. Fluid Phase Equilibria, 169, 249–276. DOI: `10.1016/S0378-3812(00)00310-1`.
13. Quiñones-Cisneros, S. E., Zéberg-Mikkelsen, C. K., Stenby, E. H. (2001). *The friction theory for viscosity modeling: extension to crude oil systems*. Chemical Engineering Science, 56, 7007–7015. DOI: `10.1016/S0009-2509(01)00335-9`.
14. Mylona, S. K., et al. (2014). *Reference Correlations for the Density and Viscosity of Squalane from 273 to 473 K at Pressures to 200 MPa*. Journal of Physical and Chemical Reference Data, 43, 013104. DOI: `10.1063/1.4863984`.
15. IAPWS (2018). *Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use*, R6-95(2018).

网络复核入口（访问日期：2026-09-16）：[Søreide–Whitson](https://doi.org/10.1016/0378-3812(92)85105-H)、[nC4–H2O 高温高压数据](https://doi.org/10.1016/j.supflu.2003.09.002)、[Teratani 重油/水 PR 参数](https://doi.org/10.1627/jpi.60.26)、[CPA 水–n-烷烃参数](https://doi.org/10.1016/j.fluid.2007.05.023)、[IAPWS-95](https://www.iapws.org/relguide/IAPWS-95.html)、[IAPWS 2008 黏度](https://www.iapws.org/relguide/viscosity.html)。
