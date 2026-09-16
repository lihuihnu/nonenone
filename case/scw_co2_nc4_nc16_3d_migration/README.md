# 三维超临界水–CO2–nC4–nC16 运移算例

本算例用于定性研究超临界水共注对轻/重烷烃运移、气体上浮、跨层窜流和重油
采出的影响。它是可复现的等温数值实验，不是某个油藏的历史拟合。

## 物理设计

| 项目 | 设置 |
|---|---:|
| 网格 | 24×16×6，共 2304 单元 |
| 尺寸 | 240×160×36 m |
| 温度 / 初始压力 | 673.15 K / 25 MPa |
| 水临界点 | 647.096 K / 22.064 MPa |
| 组分 | 独立守恒 H2O；PR 子系统 CO2、nC4、nC16 |
| 初始饱和度 O/G/W | 0.85 / 0 / 0.15 |
| 初始油相摩尔组成 CO2/nC4/nC16 | 0.02 / 0.23 / 0.75 |
| 总注入速度 | 0.025 nominal PV/year |
| 默认注入原位体积分数 SCW/CO2 | 0.75 / 0.25 |
| 生产控制 | 与注入等量的原位总体积流量 |
| 终点 | 10 年，约 0.25 nominal PV |

超临界水密度由生产代码的 IAPWS-IF97 区域分派器计算，黏度采用
McBride-Wright 纯水极限。`Water` 只是独立守恒/流动槽位名称，不表示把该流体
当作常温液态水。CO2、nC4、nC16 使用 Peng-Robinson 相平衡；nC4 代表轻烷，
nC16 代表重烷。

之所以选择“独立 SCW + 三组分 PR”而不是把 H2O 强制放入三相液-汽 Flash，是
为了在 673.15 K 的水临界点以上保持清晰的模型适用边界，并直接研究 SCW 的
密度、黏度、相对渗透率和重力项如何改变运移。当前模型不包含能量方程、热前缘、
水热裂解、焦化、改质反应或 nC16 黏度随组成变化，因此结果应解释为等温流动和
相行为效应，而不是热采采收率预测。

## 三维迁移机制

- 注入井位于西南侧下部两层，生产井位于东北侧上部三层；两井采用等量原位
  总体积率控制，避免 BHP 井在低流量比较中发生相反向窜流。
- 一条斜向高渗通道连接两井附近区域。
- `k=2` 是低渗隔层，仅在模型中央保留有限导流窗口。
- 水平/垂向各向异性、SCW/气/油密度差和上下错层完井共同形成三维迁移。

重点观察：SCW 是否延缓或促进 CO2 突破、是否改变气相质心高度、是否强化隔层
窗口处的跨层流动，以及 nC4 与 nC16 的采出差异。

## 构建与运行

```bash
make case CASE=scw_co2_nc4_nc16_3d_migration -j

# 主工况：75% SCW + 25% CO2，原位总体积流量固定
make run CASE=scw_co2_nc4_nc16_3d_migration NP=2 \
  RESULT_DIR=./results/scw-75 RUN_ARGS='-scw_fraction 0.75'

# 对照：相同原位总体积流量的纯 CO2 注入
make run CASE=scw_co2_nc4_nc16_3d_migration NP=2 \
  RESULT_DIR=./results/co2-only RUN_ARGS='-scw_fraction 0.0'
```

`-scw_fraction` 必须位于 `[0,1]`，表示注入流中 SCW 的原位相体积分数；剩余
部分为纯 CO2 气相。对照保持总原位体积流量、井位、BHP、网格和模拟时间不变。
它回答的是“以 SCW 替代部分 CO2 作为驱替介质后运移如何变化”，不能单独分离
热效应或化学反应效应。

快速预检：

```bash
case/scw_co2_nc4_nc16_3d_migration/bin/scw_co2_nc4_nc16_3d_migration \
  -numSteps 1 -dt 1e-4 -adaptive_dt false -scw_fraction 0.75 \
  -result_dir /tmp/scw-3d-smoke
```

## 对比指标

两个工况完成后运行：

```bash
python3 case/scw_co2_nc4_nc16_3d_migration/compare_scw_effect.py \
  --scw ./results/scw-75 --control ./results/co2-only \
  --output ./results/scw_effect_summary.csv
```

脚本比较最终压力跨度、平均 O/G/W 饱和度、气相覆盖范围、气相三维质心、上半部
气/水占比，以及 nC4、nC16 累计采出与初始库存之比。解释结果时还应检查：

- `component_mass_balance.csv`：逐组分守恒误差；
- `well_history.csv`：井底压力、产液/产气和控制切换；
- `solution_final.csv`：压力、饱和度和油气组成空间分布；
- `reservoir_diagnostics.csv`：全局流动与收敛诊断。

正式结论至少应报告 SCW/CO2 两个工况均通过质量守恒门槛，并对注入组成不同造成
的 CO2 质量输入差异作归一化说明。
