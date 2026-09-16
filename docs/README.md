# MPMC_SCW 文档导航

主文档只描述当前设计，历史过程由 Git 追溯。

| 文档 | 什么时候读 |
|---|---|
| [INTEGRATION_DECISIONS.md](INTEGRATION_DECISIONS.md) | 了解 W3/G8K 来源、冲突取舍与锁定基线 |
| [ARCHITECTURE.md](ARCHITECTURE.md) | 查模块边界、依赖关系与代码入口 |
| [MODEL.md](MODEL.md) | 理解未知量、方程、流动、井、时间步和输出 |
| [THERMODYNAMICS.md](THERMODYNAMICS.md) | 修改或验证 PR/SW/CPA、稳定性、Flash 与相态切换 |
| [HPC_RUN.md](HPC_RUN.md) | 配置 WSL、本地 Linux 或超算，处理 MPI/PETSc 问题 |
| [DEVELOPMENT.md](DEVELOPMENT.md) | 查看开发流程、测试矩阵、Git 与扩展位置 |
| [CONTINUE_DEVELOPMENT.md](CONTINUE_DEVELOPMENT.md) | 从交接包恢复并完成首次门禁 |
| [H2O_CO2_NC10_CASE_PARAMETERS.md](H2O_CO2_NC10_CASE_PARAMETERS.md) | 查询四模型二维 benchmark 的完整物理、井控、EOS 与数值参数台账 |
| [H2O_CO2_NC10_0P1_NUMERICAL_EXPERIMENT.md](H2O_CO2_NC10_0P1_NUMERICAL_EXPERIMENT.md) | 阅读四模型 0.1 PVI 论文式数值实验报告 |
| [SCW_KEROGEN_FLOW_EXPERIMENT.md](SCW_KEROGEN_FLOW_EXPERIMENT.md) | 超临界水–干酪根裂解产物单重质及轻/中/重拟组分驱替设计 |
| [SCW_KEROGEN_CURRENT_CASE_REPORT.md](SCW_KEROGEN_CURRENT_CASE_REPORT.md) | 旧有限初始水相参数集的历史模型、结果与优化报告 |
| [SCW_KEROGEN_2D_CALIBRATION_DESIGN.md](SCW_KEROGEN_2D_CALIBRATION_DESIGN.md) | 60×20×1 二维扩展、squalane/nC4/nC10 标定边界、黏度计划、三 EOS 差异带与优化 DOE |
| [SCW_BINARY_ZERO_FLOW_CALIBRATION.md](SCW_BINARY_ZERO_FLOW_CALIBRATION.md) | 水–角鲨烷零维 PVT、温度相关 BIP、纯组分黏度回归结果及参数接受/拒绝判据 |
| [SCW_KEROGEN_PHYSICAL_BASELINE_20260916.md](SCW_KEROGEN_PHYSICAL_BASELINE_20260916.md) | 单相起步、IAPWS 水黏度、分体系 BIP、组分见水判据及六个二维烟雾测试的当前物理基线 |
| [CHANGELOG.md](CHANGELOG.md) | 查看少量关键架构与能力里程碑 |

局部 API 优先阅读各核心目录的 `MODULE.md`。算例、测试和离线工具分别见 [../case/README.md](../case/README.md)、[../test/README.md](../test/README.md) 和 [../tools/README.md](../tools/README.md)。

巨型总手册、单次审计过程、agent 工作日志和临时修复报告不随生产仓库维护；需要时使用 `git log`、`git show` 和 tag 追溯。
