# MPMC_SCW

MPMC_SCW 是一个基于 C++17、PETSc 与 MPI 的全隐式多相多组分储层模拟器。本仓库统一了 W3 的 H2O-CO2-nC10 三相物理与 G8K 的网格、Natural 装配性能开发线，只保留一套生产实现。

当前基线包括：

- PR、Søreide-Whitson 和显式 opt-in 的 PR-CPA 热力学；
- P-T-z 稳定性分析、两/三相 Flash 与带滞回的相态切换；
- IAPWS+Garcia 富水相密度与 McBride-Wright 水相黏度；
- StructuredGridCore 与 root-only 分布式 CpGrid/GRDECL 读取；
- BHP/rate 井控、reservoir-total-rate 控制与自适应时间步；
- Natural 蓄积、面通量和 Jacobian 拓扑缓存及 PETSc 批量写入；
- 质量守恒、井控、输出与失败诊断。

生成图、运行结果和构建产物不进入 Git；可复现的绘图脚本与参考输入保留。

## 快速开始

先读 [START_HERE.md](START_HERE.md)，再为当前机器创建一个不入库的配置：

```bash
# 已验证的当前 WSL
cp config/wsl.example.mk config/hpc.local.mk

# 通用 Linux 工作站
# cp config/local.example.mk config/hpc.local.mk

make print-config
make doctor-local
make audit
make unit
```

目标超算默认值位于 `config/hpc.mk`；账号相关路径和 module 可在 `config/hpc.local.mk` 覆盖。详见 [config/README.md](config/README.md) 与 [docs/HPC_RUN.md](docs/HPC_RUN.md)。

## H2O-CO2-nC10 二维 benchmark

60x20x1 规则网格包含四种严格对齐的配置：

| 配置 | 算例 / 选择器 |
|---|---|
| Traditional | `h2o_co2_nc10_2d_traditional` |
| New-PR | `h2o_co2_nc10_2d_benchmark`, `EOS=pr` |
| New-SW | `h2o_co2_nc10_2d_benchmark`, `EOS=sw` |
| New-CPA | `h2o_co2_nc10_2d_benchmark`, `EOS=cpa` |

编译并执行一个短时步：

```bash
make case CASE=h2o_co2_nc10_2d_traditional -j
make case CASE=h2o_co2_nc10_2d_benchmark -j
make run CASE=h2o_co2_nc10_2d_benchmark NP=1 EOS=pr \
  RESULT_DIR=./results/new-pr RUN_ARGS='-numSteps 1 -dt 0.01 -adaptive_dt false'
```

## 常用命令

```bash
make help
make audit
make unit
make cases -j
make case CASE=<name> -j
make distributed-mesh NP=2
make integration NP=2 MESH_DIR=/path/to/DQ_data
make package
```

## 仓库结构

```text
ad/                    自动微分
common/                单位、数学和公共 PETSc I/O
grid/                  StructuredGrid、CpGrid、GRDECL 与分区
indices/               主变量和方程布局
models/                Natural、物性、EOS/Flash 与井
AdaptiveTimeStepper/   自适应时间步和井控循环
output/                CSV、井、守恒与 VTK 输出
case/                  正式算例与 benchmark
test/                  单元、编译和 PETSc/MPI 门禁
tools/                 离线热力学与标定工具
config/                默认配置与机器模板
docs/                  架构、物理、运行和开发文档
scripts/               审计、环境和交接工具
```

后续 agent 必须遵循 [AGENTS.md](AGENTS.md)。架构与物理语义分别见 [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)、[docs/MODEL.md](docs/MODEL.md) 和 [docs/THERMODYNAMICS.md](docs/THERMODYNAMICS.md)。
