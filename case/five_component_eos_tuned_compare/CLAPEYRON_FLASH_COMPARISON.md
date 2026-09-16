# 五组分初始闪蒸：Our 与 Clapeyron.jl 对比

## 结论

在完全相同的温度、压力、总组成和 EOS 参数下，Clapeyron.jl 0.6.27
复算得到的 PR78 与 sCPA 三相平衡均与 Our 生产程序一致。

| EOS | Clapeyron 检查方式 | 最大相分率差 | 最大饱和度差 | 最大相组成差 | 判定 |
|---|---|---:|---:|---:|---|
| PR78 | 以 Our 三相组成为初值的多相 Newton | `4.03e-6` | `3.09e-5` | `2.62e-5` | 一致 |
| sCPA | 以 Our 三相组成为初值的多相 Newton | `3.33e-16` | `5.13e-16` | `1.11e-16` | 双精度舍入误差内一致 |
| sCPA | 不使用 Our 相组成的 full-TPD 搜相 | `6.41e-9` | `2.62e-8` | `3.31e-9` | 独立搜相一致 |
| PR78 | 不使用 Our 相组成的 full-TPD 搜相 | — | — | — | Clapeyron 未返回有效三相状态 |

因此，本次交叉检查支持以下两点：

1. Our 的 PR78 和 sCPA 三相平衡方程、相组成、相分率及相体积计算与
   Clapeyron 的同参数实现相符；
2. sCPA 还通过了不依赖 Our 相组成初值的 Clapeyron full-TPD 搜相检查。

PR78 的 full-TPD 失败只说明 Clapeyron 通用自动初始化没有在该强非对称
水–烃体系中找到三相分支，不能据此否定已经收敛的三相平衡解。

## 对比条件

| 项目 | 设置 |
|---|---|
| 组分顺序 | H2O / CO2 / CH4 / C2H6 / nC4H10 |
| 总组成 | `z = [0.30, 0.10, 0.15, 0.15, 0.30]` |
| 压力 | `6.0 MPa` |
| 温度 | `305 K` |
| 目标相态 | oil + gas + water 三相 |
| Our 输入 | 重新运行生产程序得到的 `phase_state_step_0.csv` |
| 外部软件 | Clapeyron.jl `0.6.27`，Julia `1.12.7` |

没有为 Clapeyron 重新拟合参数。脚本从当前算例的 `case_config.hpp` 固化
同一套临界性质、偏心因子、纯组分 CPA 参数、BIP、4C-water、B2-CO2
以及 H2O–CO2 交叉缔合参数。

模型对应关系如下：

- Our PR 对应 Clapeyron `PR78`；
- Our CPA 对应 Clapeyron `CPA(...; radial_dist=:KG)`，即当前算例使用的
  SRK-CPA/Kontogeorgis 径向分布函数分支；
- 当前 Clapeyron 0.6.27 源码中没有 Søreide–Whitson 原生模型，因此 SW
  未纳入严格同模型对比。用普通 PR 加经验修正替代 SW 会改变模型定义，
  不能作为等价交叉验证。

## 三相结果

### PR78

| 程序 | oil 摩尔分率 | gas 摩尔分率 | water 摩尔分率 |
|---|---:|---:|---:|
| Our | 0.64888038 | 0.05142131 | 0.29969831 |
| Clapeyron | 0.64887635 | 0.05142511 | 0.29969854 |

| 程序 | oil 饱和度 | gas 饱和度 | water 饱和度 |
|---|---:|---:|---:|
| Our | 0.71754191 | 0.19959733 | 0.08286076 |
| Clapeyron | 0.71751115 | 0.19962827 | 0.08286058 |

两套结果对相态和各相占比的判断相同。最大差异出现在气相饱和度，绝对值
约 `3.09e-5`，不足 0.004 个百分点。该量级可由非线性终止阈值、相根选择
细节和物性常数末位差异解释，不影响流动初态的工程判读。

### sCPA

| 程序 | oil 摩尔分率 | gas 摩尔分率 | water 摩尔分率 |
|---|---:|---:|---:|
| Our | 0.67723622 | 0.02227628 | 0.30048750 |
| Clapeyron | 0.67723622 | 0.02227628 | 0.30048750 |

| 程序 | oil 饱和度 | gas 饱和度 | water 饱和度 |
|---|---:|---:|---:|
| Our | 0.82385781 | 0.09889628 | 0.07724591 |
| Clapeyron | 0.82385781 | 0.09889628 | 0.07724591 |

定向多相 Newton 的全部相分率、饱和度、压缩因子和 15 个相组成值在
`1e-16` 量级内相同。无相组成初值的 full-TPD 搜索也独立找到同一三相根，
其最大饱和度差仅为 `2.62e-8`。

Clapeyron 对当前非对称 H2O–CO2 B2 交叉缔合拓扑必须启用
`AssocOptions(implicit_ad=true)`。不开启时，缔合位点分数本身可以求得，
但压力导数会出现 `NaN`；开启后，按 Our 三相状态回算的三相压力均为
`6.000000 MPa`，随后 sCPA 闪蒸正常收敛。

## 相组成与压缩因子

逐相、逐组分的完整值保存在
[`clapeyron_comparison/phase_compositions.csv`](clapeyron_comparison/phase_compositions.csv)，
相分率、饱和度和压缩因子保存在
[`clapeyron_comparison/initial_flash_summary.csv`](clapeyron_comparison/initial_flash_summary.csv)。

最大绝对差汇总为：

| EOS | `max |Δbeta|` | `max |ΔS|` | `max |ΔZ|` | `max |Δx_phase|` |
|---|---:|---:|---:|---:|
| PR78 | 4.0308e-6 | 3.0939e-5 | 1.9688e-5 | 2.6231e-5 |
| sCPA | 3.3307e-16 | 5.1348e-16 | 1.1102e-16 | 1.1102e-16 |

这里的饱和度不是 `tp_flash` 直接输出，而是两端各自根据相摩尔分率和相摩尔
体积计算：`S_p = beta_p V_p / sum(beta_q V_q)`。Clapeyron 中 oil/water
使用液相体积根，gas 使用气相体积根。

## 运行时间如何理解

本次独立 Julia 进程中，定向 PR78 和 sCPA 闪蒸的墙钟时间分别为
`41.3 s` 与 `202.3 s`。这些数字包含 Clapeyron/Julia 首次方法编译，不能与
已编译 C++ 生产程序直接比较，也不能当作稳态单次闪蒸耗时。CSV 中 Our 的
时间写为 `NaN`，表示本次没有在同一计时边界内测量，而不是零耗时。

本实验的主要评价量是同一热力学状态的数值一致性，不是跨语言冷启动性能。
若要比较求解性能，应在两个程序中分别预热后，对同一批 P–T–z 状态重复
闪蒸并统一统计中位数、P95 和失败率。

## 可复现方法

在 WSL 中执行定向三相对比：

```bash
export JULIA_DEPOT_PATH=/home/mpmc/.local/share/julia-depot
cd case/five_component_eos_tuned_compare
/home/mpmc/.local/bin/julia \
  --project=/home/mpmc/.local/share/clapeyron-mpmc \
  clapeyron_flash_compare.jl
```

增加不使用 Our 相组成初值的 full-TPD 检查：

```bash
CLAPEYRON_FULL_TPD=1 /home/mpmc/.local/bin/julia \
  --project=/home/mpmc/.local/share/clapeyron-mpmc \
  clapeyron_flash_compare.jl
```

脚本会重新生成三个主 CSV；启用 full-TPD 时还会生成
`full_tpd_check.csv`。脚本只读取本算例结果，不修改生产 EOS 或闪蒸代码。

## 适用边界

- 这是当前单一 `P–T–z` 初始状态的交叉实现检查，不等同于整个相包络的验证；
- 定向 Newton 证明两端在同一三相根上满足一致的平衡方程，但其相态初值来自
  Our；full-TPD 结果单独列出，避免把局部复算误称为全局独立搜相；
- SW 因缺少 Clapeyron 原生等价模型而未强行比较；
- PR 的 full-TPD 自动搜相失败应作为外部软件初始化能力限制记录，而不是
  擅自改用二相结果。

## 软件依据

- Clapeyron.jl：<https://github.com/ClapeyronThermo/Clapeyron.jl>
- `tp_flash` 文档：<https://clapeyronthermo.github.io/Clapeyron.jl/dev/properties/flash/>
- D.-Y. Peng, D. B. Robinson. A New Two-Constant Equation of State.
  *Industrial & Engineering Chemistry Fundamentals*, 1976, 15(1): 59–64.
- G. M. Kontogeorgis et al. An Equation of State for Associating Fluids.
  *Industrial & Engineering Chemistry Research*, 1996, 35(11): 4310–4320.
