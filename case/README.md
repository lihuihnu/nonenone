# case 模块

`case/` 是正式储层算例的唯一入口，只负责组装生产模块，不复制 EOS、Flash、网格、井或求解器实现。

## 目录约定

```text
case/<name>/
├── <name>.cpp          主程序装配
├── case_config.hpp     模型、流体、网格、初值、时间与输出
├── well_config.hpp     井与 schedule
├── case_fluid.hpp      可选：算例专属闭包或表格
├── README.md           可选：物理来源与运行说明
└── scripts/            可选：可复现后处理源码
```

运行生成的 `bin/`、`logs/`、`results/` 和 `run.sh` 不提交。

## 当前算例

| 算例 | 主要用途 |
|---|---|
| `3p6c_original` | 旧三相六组分基线与兼容路径 |
| `DQcase` | GRDECL/CpGrid、真实网格与并行路径 |
| `3p4c_pr_reservoir` | 四组分 fully-compositional PR |
| `ma2021_5c_three_phase_reservoir` | 五组分三相文献型 benchmark |
| `panfili2025_case3_fullphysics` | 文献条件下的 full-physics 对照 |
| `three_eos_3d_compare` | 统一三维条件下 PR/SW/CPA 对比 |
| `five_component_eos_tuned_compare` | 标定五组分三 EOS 对比 |
| `five_component_eos_tuned_legacy_compare` | traditional independent-water + PR 对照 |
| `h2o_co2_nc10_2d_traditional` | 60x20x1 Traditional 基线 |
| `h2o_co2_nc10_2d_benchmark` | 相同控制量下 New-PR/New-SW/New-CPA 对比 |
| `h2o_co2_nc10_2d_scw_traditional` | 60x20x1 Traditional 超临界水对照 |
| `h2o_co2_nc10_2d_scw_benchmark` | 60x20x1 New-PR/New-SW/New-CPA 超临界水对比 |
| `h2o_co2_binary_2d_conventional_traditional` | 常规温压 H2O-CO2 Traditional 二元基线 |
| `h2o_co2_binary_2d_conventional_new` | 常规温压 H2O-CO2 New-PR/New-SW/New-CPA 二元对比 |
| `h2o_co2_binary_2d_scw_traditional` | 超临界水条件 H2O-CO2 Traditional 二元基线 |
| `h2o_co2_binary_2d_scw_new` | 超临界水条件 H2O-CO2 New-PR/New-SW/New-CPA 二元对比 |
| `sun2024_exp12_scw_co2_nc16_1d` | Exp.12 派生独立超临界水 + PR(CO2-nC16) 一维驱替 |
| `sun2024_scw_co2_nc16_factorial_1d` | Sun 2024 Exp.4/6/10/12 2×2复现矩阵及等总流量对照 |
| `h2o_co2_bsb_lumped_3d_lab` | 复用 60×20×1 nC10 网格的 BSB 代表组分 CO2/SCW 五点驱替对照 |
| `scw_co2_nc4_nc16_3d_migration` | 异质规则网格中 SCW/CO2 对 nC4/nC16 三维运移影响 |
| `scw_kerogen_squalane_1d` | H2O–squalane 单重质拟组分严格 SCW 一维驱替 |
| `scw_kerogen_lmh_1d` | H2O + 轻/中/重裂解产物拟组分严格 SCW 一维驱替 |
| `scw_kerogen_squalane_2d` | 60×20×1 H2O–squalane PR/SW/CPA 二维模型差异带 |
| `scw_kerogen_lmh_2d` | 60×20×1 H2O + 轻/中/重组分 PR/SW/CPA 二维模型差异带 |

真实参数以各目录的 `case_config.hpp` 和 `well_config.hpp` 为准。
二元算例的统一运行与绘图清单见
[`h2o_co2_binary_unified_compare/README.md`](h2o_co2_binary_unified_compare/README.md)。

## 编译与运行

```bash
make case CASE=h2o_co2_nc10_2d_benchmark -j
make run CASE=h2o_co2_nc10_2d_benchmark NP=2 EOS=sw RESULT_DIR=./results/sw
```

也可在本目录运行：

```bash
make list
make h2o_co2_nc10_2d_traditional -j
make h2o_co2_nc10_2d_benchmark -j
make prepare CASE=h2o_co2_nc10_2d_benchmark
```

`prepare` 按当前 machine profile 生成带管理标记的 `run.sh`。机器路径、launcher 和 module 只放在 `config/hpc.local.mk`；额外 PETSc 参数通过命令行传入。

## 新增算例

从最接近的算例复制目录并重命名主程序，然后只修改装配、配置和必要的专属闭包。`case/Makefile` 会自动发现含 `case_config.hpp` 的目录，无需维护手写目标清单。

长期模型语义见 [../docs/MODEL.md](../docs/MODEL.md)，热力学约束见 [../docs/THERMODYNAMICS.md](../docs/THERMODYNAMICS.md)。
