# MPMC_SCW 关键里程碑

这里只记录会改变“当前项目如何理解”的重要里程碑。逐次 bug 修复、性能调优和实验过程请直接查看 Git：

```bash
git log --oneline --decorate --all
```

## 2026-09 — 超临界水相态恢复与统一对比算例

- 加固 SW 限制三相 Flash 的温度路径恢复、相角色归一化与高温稳定性回归；
- 扩展超临界水区域的 Garcia 溶解 CO2 偏摩尔体积闭合，并保留锁定的 IAPWS 水基准；
- 增加常规/超临界条件下 H2O-CO2 二元与 H2O-CO2-nC10 三元对齐算例；
- 增加 Sun 2024 一维驱替、SCW/CO2/nC4/nC16 三维运移和 Heringer 2025 三相闪蒸基准。

## 2026-08 — 可复现交接与网格接口加固

- 标准交接包改为 committed source snapshot + `repository.bundle` + manifest + SHA-256；
- 本机配置、外部网格、结果和未提交工作树不进入交接包；
- StructuredGrid/CpGrid/GRDECL 输入与几何接口继续加固；
- DQ 与多套真实 corner-point 网格进入并行验证路径。

## 2026-08 — 五组分 EOS 对比体系

- 建立 fully-compositional PR / SW / CPA 五组分对比算例；
- 增加传统 independent-water + PR 对照算例；
- 增强公开数据和独立 Flash/EOS 交叉验证。

## v67.12 — Flash/phase validation 与 PETSc hook 兼容

- PR/SW/CPA 共用 fail-only pressure continuation 基础设施；
- 标准 PETSc 使用公开 line-search callback；
- 保留历史定制 PETSc `updateState/updateSol` C/C++ ABI 兼容；
- Flash/phase validation 从 reservoir case 中进一步独立。

## v67.10–v67.11 — Runtime、case 与运行层收敛

- Natural runtime/state 责任进一步拆分；
- 正式 case 统一井 factory 与配置方式；
- 本地 Linux/WSL 与 HPC 使用统一 profile/Makefile 工作流。

## v67.5–v67.9 — CPA、active-set 与代码结构整理

- CPA calibration/holdout 与性能路径加强；
- fully-compositional phase active-set、事务式状态和 failure channel 固化；
- 大型 Natural/Flash 头文件逐步拆分；
- 测试和工具入口统一。

## v59–v67.3 — 三 EOS 三相互溶生产基础

- O/G/W 固定 public slot 与内部 thermodynamic role 解耦；
- SW 相边界恢复、压力延续与 restricted/full Flash 路径形成；
- PR/SW/CPA 三种 backend 纳入统一三相 Flash 框架；
- Natural 三相多组分主模型成为后续生产基线。
