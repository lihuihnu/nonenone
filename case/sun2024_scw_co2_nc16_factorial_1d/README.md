# Sun-2024超临界水–CO2砂管可复现实验矩阵

本目录把 Sun et al. (2024), DOI `10.3390/app14093588` 中同为400 °C的
Exp.4、6、10、12组织为一个预注册的2×2数值实验，并增加两个等储层原位体积流量
对照。它只使用论文物理输入和生产代码，不读取或依赖任何已有数值结果。

## 科学设计

主因子为压力（23/24 MPa）与CO2加入量（0/2 mL/min）。原实验加入CO2时，
总报告流量由10增至12 mL/min，因此R10-R04或R12-R06不能单独解释为组成效应。
C23/C24使用纯SCW，并在选定流量基准后自动匹配R10/R12的**储层总原位流量**，
专门估计流量增加造成的差异。只有在`in_situ_volume`解释下，对照的等效水量才
恰好是12 mL/min；参考密度换算时由程序反算，不能继续写死为12。

| ID | 来源 | 压力 | SCW | CO2 | 角色 |
|---|---:|---:|---:|---:|---|
| R04 | Exp.4 | 23 MPa | 10 mL/min | 0 | 物理复现 |
| R06 | Exp.6 | 24 MPa | 10 mL/min | 0 | 物理复现 |
| R10 | Exp.10 | 23 MPa | 10 mL/min | 2 mL/min | 物理复现 |
| R12 | Exp.12 | 24 MPa | 10 mL/min | 2 mL/min | 物理复现 |
| C23 | 新增 | 23 MPa | 自动反算 | 0 | R10等储层流量对照 |
| C24 | 新增 | 24 MPa | 自动反算 | 0 | R12等储层流量对照 |

完整输入登记见`experiment_manifest.csv`，实验实测终点单独保存在
`physical_validation_targets.csv`。数值参数冻结后才应把预测与实测终点比较。

## 当前模型边界

- 0.48 m×0.039 m砂管，孔隙度0.39，渗透率2000 mD；
- 初始S(O/G/W)=0.935/0/0.065；
- 一维等温模型，固定673.15 K；
- H2O为独立守恒SCW相，密度使用IAPWS-IF97，黏度使用McBride-Wright；
- CO2-nC16使用PR78与`kij=0.09`；
- 不包含能量方程、热损失、反应、改质、毛管压力、扩散或有限速率传质。

因此，热利用率不是本算例的验收量；在加入能量方程和真实重油之前，也不能把
终点误差全部解释为相对渗透率误差。

## 原文核验与流量基准

原文图1给出了柱塞泵、高压中间容器、出口背压阀（BPR）和35 MPa的超临界水
发生器上限；表2给出10 mL/min水和0或2 mL/min CO2。但是论文没有报告这两个
体积率的计量温度、计量压力或密度。方法部分还同时写了“温度稳定后停止”、
“每轮2 h”，而结果图统一以0–4 PV展示；因此本算例把4 PV定义为**共同报告终点**，
不把它声称为已完全确认的实验停机规则。详细判读见`SOURCE_AUDIT.md`。

运行时必须在以下两种基准中显式选一种：

- `in_situ_volume`：直接把表2体积率解释为673.15 K、对应背压下的原位体积率；
  这是可重复的敏感性边界，不应标为严格物理复现。
- `reference_density`：把表2体积率分别乘以实验计量状态的水/CO2密度得到质量率，
  再除以本模型在673.15 K、对应压力下的相密度得到储层体积率。这是拿到实验台账
  后的首选复现路径，且必须提供`-water_reference_density_kg_m3`和
  `-co2_reference_density_kg_m3`。

每次运行会在`experiment_metadata.csv`中记录报告体积率、参考密度、模型储层密度、
换算后的两相质量率与原位体积率，避免同名“mL/min”被混用。

## 构建与运行

```bash
make case CASE=sun2024_scw_co2_nc16_factorial_1d -j

make run CASE=sun2024_scw_co2_nc16_factorial_1d NP=1 \
  RESULT_DIR=./results/sun2024-factorial/R12 \
  RUN_ARGS='-experiment_id R12 -rate_basis in_situ_volume \
            -nx 48 -target_pv 4 -dt_pv 0.02'
```

拿到实验计量状态密度后，正式复现命令为：

```bash
make run CASE=sun2024_scw_co2_nc16_factorial_1d NP=1 \
  RESULT_DIR=./results/sun2024-factorial/reference/R12 \
  RUN_ARGS='-experiment_id R12 -rate_basis reference_density \
            -water_reference_density_kg_m3 WATER_DENSITY \
            -co2_reference_density_kg_m3 CO2_DENSITY \
            -nx 48 -target_pv 4 -dt_pv 0.02'
```

其中`WATER_DENSITY`与`CO2_DENSITY`必须来自实验泵/中间容器的实际计量状态，
不能用储层密度或随意选取的标准态密度替代。

这里禁止混用通用的`-numSteps/-dt`，终点和输出间隔统一使用PVI表达。快速预检：

```bash
case/sun2024_scw_co2_nc16_factorial_1d/bin/sun2024_scw_co2_nc16_factorial_1d \
  -experiment_id R12 -rate_basis in_situ_volume \
  -nx 48 -target_pv 1e-4 -dt_pv 1e-4 -adaptive_dt false \
  -result_dir /tmp/sun2024-factorial-smoke
```

`scripts/run_matrix.sh`可生成六工况主矩阵以及R12的网格/时间步收敛组。脚本要求
显式环境变量，避免错误基准下的批量运行：

```bash
RATE_BASIS=in_situ_volume \
  bash case/sun2024_scw_co2_nc16_factorial_1d/scripts/run_matrix.sh smoke

RATE_BASIS=reference_density \
WATER_REFERENCE_DENSITY_KG_M3=MEASURED_WATER_DENSITY \
CO2_REFERENCE_DENSITY_KG_M3=MEASURED_CO2_DENSITY \
  bash case/sun2024_scw_co2_nc16_factorial_1d/scripts/run_matrix.sh formal
```

正式矩阵默认采用已通过收敛检查的`nx=96`与`dt_pv=0.01`；可用
`FORMAL_NX`和`FORMAL_DT_PV`显式覆盖，覆盖值会写入每个结果目录的运行环境记录。

## 数值验收

先完成R12的收敛检查，再执行全部主矩阵：

```bash
RATE_BASIS=in_situ_volume \
  bash case/sun2024_scw_co2_nc16_factorial_1d/scripts/run_matrix.sh grid
RATE_BASIS=in_situ_volume \
  bash case/sun2024_scw_co2_nc16_factorial_1d/scripts/run_matrix.sh time

python3 case/sun2024_scw_co2_nc16_factorial_1d/scripts/check_convergence.py \
  --coarse ./results/sun2024-factorial/grid/R12_nx48 \
  --medium ./results/sun2024-factorial/grid/R12_nx96 \
  --fine ./results/sun2024-factorial/grid/R12_nx192 \
  --output ./results/sun2024-factorial/grid/convergence_summary.csv
```

门槛为：以`max(初始库存, 当前库存, 累计注入, 累计产出)`归一化的最大组分
守恒闭合误差`1e-6`，中/细网格nC16终点采收率差不超过
0.005，含CO2工况的突破时刻差不超过0.02 PV。CO2突破定义为离散输出区间内
“出口CO2组分质量率/入口CO2质量率”首次达到1%；若两组均未在4 PV内突破，
突破时刻差记为0，但“未突破”本身必须作为物理验证差异报告。脚本读取每次新运行生成的
`experiment_metadata.csv`，不会硬编码已有数值结果。

该砂管各组分库存远小于1 kg，不能直接使用输出CSV中按全局单位尺度保护的
`relative_error`列作为验收量；检查脚本会由`balance_error_kg`重新计算无量纲闭合。
算例还默认设置`-snes_atol 1e-12`，防止约化残差略低于PETSc默认绝对阈值时
零迭代接受时间步、造成状态冻结。用户显式传入的`-snes_atol`仍优先。

## 科研图与算例报告

在六工况正式矩阵、网格收敛组和时间步收敛组完成后，可由指定的新结果目录生成
科研图、机器可读汇总和中文技术报告：

```powershell
python case/sun2024_scw_co2_nc16_factorial_1d/scripts/generate_report_artifacts.py `
  --formal-root tmp/sun2024_factorial_report_v1/formal `
  --case-root case/sun2024_scw_co2_nc16_factorial_1d `
  --grid-summary tmp/sun2024_factorial_convergence_v4/grid/convergence_summary.csv `
  --time-summary tmp/sun2024_factorial_convergence_v4/time/convergence_summary.csv `
  --output-root output/sun2024_factorial_report
```

`generate_report_artifacts.py`只读取命令行明确给出的结果目录，不搜索已有
`results`目录；所有曲线保持原始离散点，不进行平滑。图件数据和变换记录写入
`output/sun2024_factorial_report/data/`。

## 解释顺序

1. R06-R04给出压力变化下的总差异；
2. R10-C23与R12-C24给出同一声明基准下、等总原位体积率的注入组成差异；
3. C23-R04与C24-R06给出纯SCW的流量效应；
4. 参数锁定后再与物理实验终点比较；
5. 本算例只验证独立SCW流动，不能替代PPT中H2O跨三相分配的BSB机理算例。
