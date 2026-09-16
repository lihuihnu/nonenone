# MPMC_SCW 物理与求解模型

本文只描述储层求解主路径中长期有效的模型语义。热力学公式与 Flash 细节见 [THERMODYNAMICS.md](THERMODYNAMICS.md)，软件依赖关系见 [ARCHITECTURE.md](ARCHITECTURE.md)。

## 1. 两种相行为模式

编译期配置通过 `CompositionalModelConfig` 选择相行为。

### LegacyOilGasWithIndependentWater

- 油/气使用 compositional EOS；
- 水作为独立水相；
- 可叠加旧式水相 CO2 溶解、吸附、Land trapping 等可选机制；
- 主要用于兼容既有算例。

### FullyCompositionalThreePhase

- H2O 与 CO2、CH4、烃组分统一进入组分集合；
- 三个公共相槽固定为 O / G / W；
- 所有组分可在三相之间重新分配；
- 相出现/消失不改变全局未知量维数。

## 2. Fully-compositional 未知量

对 `N` 个组分，单元主变量保持固定布局：

```text
p
x_o[0 ... N-2]
x_g[0 ... N-2]
x_w[0 ... N-2]
S_o
S_g
S_w
[p_bhp]   # 仅当模型包含井未知量
```

每相最后一个摩尔分数由归一化闭合得到，因此每相只显式保存 `N-1` 个组成未知量。

相消失时不删除未知量，而由 phase active-set 与 inactive-equation replacement 保持系统维数固定。

## 3. 方程布局

Fully-compositional 主系统由以下方程组成：

1. `N` 个组分质量守恒；
2. `N` 个 O-G fugacity equilibrium；
3. `N` 个 O-W fugacity equilibrium；
4. 1 个 saturation closure；
5. 可选 1 个 well-control equation。

active phase set 改变时，对应平衡方程会被相态约束替换，但未知量位置不变。

## 4. 组分质量守恒

每个组分满足离散形式：

$$
\frac{M_i^{n+1}-M_i^n}{\Delta t}
+ \sum_f F_{i,f}
- Q_i = 0.
$$

其中：

- $M_i$：单元内组分库存；
- $F_{i,f}$：穿过面的组分质量通量；
- $Q_i$：井或其它源汇项。

Fully-compositional 模式中 H2O 就是普通 conserved component，不再另设独立“水质量方程”。

## 5. Flash 到饱和度

P–T–z Flash 给出每相组成与相摩尔分率 `beta`。EOS 再给出各相摩尔密度，因此相体积与饱和度由：

```text
P,T,z
 -> phase split / composition
 -> EOS molar density
 -> phase molar volume
 -> phase volume
 -> saturation
```

计算得到。

`beta` 是摩尔分率，`S` 是孔隙体积分率，二者只有在特殊条件下才近似相等。

Volume translation 若启用，只进入相体积/密度相关量，不改变 fugacity、stability 或 phase split。

## 6. 达西流动

各相使用 Darcy/TPFA 形式的面流率。相势差包含压力与重力项，面 mobility 采用确定性的上游规则。

组分通量由各相流率乘相内组分含量求和：

$$
F_i = \sum_{\alpha} q_\alpha\,\rho_{m,\alpha}\,x_{\alpha i},
$$

实际实现还包含单位和质量/摩尔量之间的统一换算。

相对渗透率、黏度和密度共同决定 phase mobility。

## 7. 自动微分与 Jacobian

Natural 层使用前向 AD 传播局部残差的一阶导数。PETSc runtime 负责：

- owned/ghost 状态同步；
- 单元与邻居 block Jacobian 装配；
- SNES nonlinear solve；
- KSP/PC 线性求解；
- failure/domain-error 与时间步控制的接口。

物理层不直接承担 MPI 生命周期或全局矩阵管理。

## 8. 相态切换

公共相槽共有 7 种非空组合：

```text
O
G
W
O+G
O+W
G+W
O+G+W
```

典型更新逻辑：

```text
当前 primary + phase state
 -> 检查 saturation / active-set
 -> 候选相消失或缺失相稳定性测试
 -> restricted/full flash
 -> canonical O/G/W state
 -> 事务式提交
```

失败的候选相态不会部分写入正式状态；时间步接受后才成为新 history。

## 9. 井模型

井模块定义井几何、完井、控制和 schedule；Natural 只消费统一井接口。

当前主要控制：

- BHP；
- TotalRate；
- OilRate；
- GasRate；
- WaterRate。

符号约定：

```text
注入 > 0
采出 < 0
```

控制残差：

```text
BHP:  p_bhp - p_target = 0
RATE: q_control - q_target = 0
```

井控切换发生在同一个物理 timestep 内，并触发该 timestep 的重新求解，而不是等待下一输出步。

## 10. AdaptiveTimeStepper

时间推进具有事务语义：

```text
history state
 -> 创建 attempt
 -> nonlinear solve
 -> well-control cycle
 -> accept 或 reject
```

接受时提交：

- primary state；
- secondary/phase state；
- history；
- 井状态；
- 时间与统计。

拒绝时 rollback，并按策略缩小 `dt`。成功且迭代较容易时可增长 `dt`，但固定输出时间不会被跨越。

## 11. 网格

### StructuredGrid

规则笛卡尔网格，适合受控 benchmark、三 EOS 对比和开发回归。

### CpGrid / GRDECL

用于角点/真实油藏网格，负责：

- COORD / ZCORN / ACTNUM；
- PORO / PERM；
- MAPAXES；
- active-cell topology；
- 几何与 transmissibility 所需数据；
- METIS partition；
- owned/ghost 和 DOF 映射。

上层 Natural 通过 backend 适配接口使用网格，不应依赖具体 parser 或内部 mesh 类型。

## 12. 输出与诊断

标准输出主要包括：

- solution/state CSV；
- phase-state CSV；
- well history；
- component mass series / mass-balance ledger；
- simulation summary；
- 可选 VTK。

输出层只读已接受状态，不参与 Newton 更新或相态切换。

## 13. 可选旧物理

Legacy 模式仍可能启用：

- 独立水相 CO2 dissolution；
- competitive adsorption；
- Land residual-gas trapping。

这些机制与 fully-compositional 三相互溶不是同一物理路径。新增研究算例必须明确自己使用哪一种模型，不允许同时以两套机制重复表达同一相间传质。

## 14. 不可破坏约束

修改求解器时优先保持：

1. 固定未知量/方程布局与 phase active-set 语义一致；
2. Flash 失败不能被伪装成收敛；
3. 状态更新必须事务化，失败 attempt 不污染 history；
4. 组分守恒和井符号约定不因 backend 改变；
5. Tools/测试复用 production 物理实现，不复制第二套 EOS/Flash；
6. 网格、井、热力学、时间推进和输出保持职责分离。
