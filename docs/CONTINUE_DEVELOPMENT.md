# MPMC_SCW 交接快速入口

本文保留为标准交接包的 `START_HERE.md` 来源，只提供最短接手路径，不再承担完整开发手册职责。

## 1. 恢复仓库

标准交接包包含 committed source snapshot、`repository.bundle`、manifest 和 SHA-256 校验。

正式继续开发优先从 bundle 恢复：

```bash
git clone repository.bundle MPMC_SCW_case
cd MPMC_SCW_case
git status --short --branch
```

`source/` 适合快速查看或临时编译；长期开发应使用 Git 仓库。

## 2. 新机器配置

```bash
cp config/local.example.mk config/hpc.local.mk
# 修改本机 PETSc / MPI / FMT / METIS 路径
make print-config
make doctor-local
```

超算环境见 [HPC_RUN.md](HPC_RUN.md)。

## 3. 首轮门禁

```bash
make audit
make unit
make cases -j
```

PETSc/MPI/真实网格配置完成后，再执行对应 integration/full 门禁。

## 4. 阅读顺序

```text
README.md
  -> docs/ARCHITECTURE.md
  -> docs/MODEL.md
  -> docs/THERMODYNAMICS.md   # 修改热力学时
  -> 对应模块 MODULE.md
```

开发、Git、测试矩阵与交接规则统一见 [DEVELOPMENT.md](DEVELOPMENT.md)。

## 5. 重要原则

- 当前状态以 Git `HEAD`、源码和实际测试结果为准；
- 不从旧压缩包拼接代码；
- 不提交 `config/hpc.local.mk`、运行结果、日志或机器私有路径；
- 不通过放松物理残差掩盖 Flash/SNES/MPI 错误；
- case/test/tools 复用 production API，不复制第二套实现；
- 历史过程通过 Git log/tag 查询，不在主文档重复保存。
