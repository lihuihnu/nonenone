# Sun et al. (2024) Exp. 12 派生超临界水–CO2–nC16 一维驱替

本算例以 Sun 等公开的一维砂管实验 Exp. 12 为几何和注采工况原型，使用
生产 Peng–Robinson 路径计算 CO2–nC16 相态，并把超临界水作为独立、可流动、
可压缩的水连续相。它是可复现的**等温派生算例**，不是原始重油热反应实验的
历史拟合。

## 公开实验输入

| 参数 | Exp. 12 |
|---|---:|
| 砂管长度 / 内径 | 0.48 m / 0.039 m |
| 孔隙度 / 渗透率 | 0.39 / 2000 mD |
| 初始油饱和度 | 0.935 |
| 温度 / 压力 | 673.15 K / 24 MPa |
| 水 / CO2 流量 | 10 / 2 mL min⁻¹ |
| 实验终点驱油效率 | 82.96% |

来源：Yan Sun et al., *Applied Sciences* 14 (2024) 3588,
[DOI 10.3390/app14093588](https://doi.org/10.3390/app14093588)。水临界点取
IAPWS 的 647.096 K、22.064 MPa。

## 模型闭合与适用边界

- Natural 当前为等温模型，全域固定在 673.15 K、24 MPa；未模拟原实验从
  50 °C 开始的升温前缘、热损失、水热裂解或改质反应。
- 真实 Bohai 重油以 nC16 代替。CO2–nC16 使用 PR78 和显式 `kij=0.09`；
  超临界水不参与油气闪蒸，作为独立守恒相处理。
- 独立水相不是“把超临界水当液态水”：其 24 MPa、673.15 K 密度由
  IAPWS‑IF97 Region 2 计算（148.5610835 kg m⁻³），黏度采用
  McBride–Wright 纯水极限。框架中的 `W` 仅是守恒/流动槽位名称。
- 论文的 48×19×19 网格被体积保持地约化为 48×1×1；轴向仍为 1 cm。
- 单入口井同时注入 5/6 体积分数的超临界水相和 1/6 体积分数的纯 CO2 相，
  总原位体积流量为 12 mL min⁻¹。
- 相对渗透率数据未公开，因此使用 `benchmark_common.hpp` 中显式记录的 Corey
  工程曲线。终点取累计注入 4 PV，也属于派生设置。
- 82.96% 只作实验背景。只有加入能量方程、真实重油组分及反应动力学后，才能
  将本算例升级为定量历史拟合。

## 构建与运行

```bash
make case CASE=sun2024_exp12_scw_co2_nc16_1d -j

make run CASE=sun2024_exp12_scw_co2_nc16_1d NP=1 \
  RESULT_DIR=./results/sun2024/new-pr
```

快速预检：

```bash
case/sun2024_exp12_scw_co2_nc16_1d/bin/sun2024_exp12_scw_co2_nc16_1d \
  -numSteps 1 -dt 1e-6 -adaptive_dt false \
  -result_dir /tmp/sun2024-pr-smoke
```

主要验证量包括 nC16 累计产出、CO2 突破、压降、O/G/W 饱和度前缘和逐组分
质量守恒。nC16 采收率可由 `component_mass_balance.csv` 计算：

```text
RF_nC16 = cumulative_produced_kg(nC16) / initial_inventory_kg(nC16)
```
