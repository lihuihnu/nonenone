# MPMC_SCW 开发与交接规则

本文只保存当前工作规则，不记录某次 agent 的过程日志或已经完成的审计报告。

## 1. 开始开发

先确认基线：

```bash
git status --short --branch
git log -5 --oneline --decorate
```

独立任务使用独立分支：

```bash
git switch -c feature/<topic>
```

已有未提交修改时先识别范围，不使用 `git reset --hard` 或批量 checkout 覆盖未知工作。

## 2. 修改原则

1. production 模块保持清晰单向依赖；
2. case 不复制 EOS、Flash、井或网格实现；
3. tools/test 尽量调用 production API；
4. 物理修改、数值恢复、机械重构、文档移动分开提交；
5. 输入错误应尽早失败；
6. 不通过放松守恒/平衡残差来伪造收敛；
7. MPI collective 使用对象自己的 communicator；
8. 本机路径、结果、日志、二进制和可重建数据不提交。

## 3. 新功能放在哪里

| 功能 | 位置 |
|---|---|
| EOS / stability / Flash | `models/include/natural/thermo/` |
| 密度、黏度、相对渗透率等闭合 | `models/include/natural/properties/` |
| accumulation / flux / source | `models/include/natural/physics/` |
| PETSc runtime / assembly | `models/include/natural/petsc/` |
| 井模型 / 控制 / schedule | `models/include/well/` |
| Structured/CpGrid/GRDECL | `grid/` |
| 时间步策略 | `AdaptiveTimeStepper/` |
| 输出和诊断 | `output/` |
| 正式算例 | `case/<name>/` |
| 回归测试 | `test/` |
| 相图、PVT、离线扫描 | `tools/` |

新增跨模块“万能 helper”前先判断它到底属于哪个责任层。

## 4. 测试矩阵

最低要求按影响范围选择，不机械跑无关测试：

| 修改范围 | 最低门禁 |
|---|---|
| 仅文档 | `python3 scripts/check_markdown_links.py` + `make audit` |
| 构建/打包脚本 | 语法检查 + 相关命令实跑 + `make audit` |
| PETSc-free production | `make unit` |
| EOS/Flash/物性 | `make unit` + 对应 tools/public-data regression + case preflight |
| StructuredGrid/PETSc runtime | unit + cases + 1/2/4 rank integration |
| CpGrid/GRDECL/partition | unit + DQ 1/2/4 rank + 至少一种其它真实网格 |
| 井/时间步/输出 | 对应 unit + Natural integration + representative case |

稳定基线可额外执行：

```bash
make unit-gcc
make unit-clang
```

测试记录至少包含：命令、编译器、MPI rank、数据集、误差/结果和未执行项。编译成功不能替代运行成功。

## 5. 文档规则

项目级长期文档只保留：

```text
README.md
docs/README.md
docs/ARCHITECTURE.md
docs/MODEL.md
docs/THERMODYNAMICS.md
docs/HPC_RUN.md
docs/DEVELOPMENT.md
docs/CHANGELOG.md
```

局部接口写就地 `MODULE.md`，算例细节写对应 case README。

不要新增：

```text
*_Vxx.md
*_FIX_REPORT.md
*_AUDIT_REPORT.md
agent_worklog.md
```

当前设计写主题文档；历史过程写 Git commit。

## 6. Commit 规则

推荐前缀：

```text
fix(scope): ...
refactor(scope): ...
test(scope): ...
docs(scope): ...
chore(scope): ...
```

提交前：

```bash
git diff --check
git status --short
git add -- <明确文件列表>
git diff --cached --check
```

避免 `git add -A` 无意带入结果、本机 profile 或别人修改。

内部 stable tag 不移动、不覆盖。

## 7. Generated data

默认不跟踪：

```text
build/
bin/
logs/
results/
reference_output/  # 若可确定性重建
临时 MATLAB/plot 输出
本地缓存
config/hpc.local.mk
```

只有不可再生且属于功能必需的小型 fixture 才进入仓库。

## 8. VS Code / IntelliSense

仓库只维护共享 C++ IntelliSense 所需配置；个人 `settings.json`、`tasks.json`、`launch.json` 不进入项目。

修改默认 PETSc/FMT/METIS/VTK include 路径或项目 include 目录后，应同步 `.vscode/c_cpp_properties.json`，并执行 `make audit`。

IntelliSense 不是构建真相；实际编译仍以 Makefile 为准。

## 9. 源码注释

注释主要解释：

- 接口契约；
- 单位与符号；
- 非显然物理/数学关系；
- 数值恢复为何存在；
- 状态与并行约束。

不要给显然 getter、循环、标准库调用或直接赋值添加噪声注释。修改算法时同步修改附近注释。

## 10. 标准交接包

生成：

```bash
make package
```

验证：

```bash
make verify-package PACKAGE=/path/to/MPMC_SCW_case_v67_12-handoff-<commit>.zip
```

包应包含：

```text
START_HERE.md
HANDOFF_MANIFEST.txt
SHA256SUMS
repository.bundle
source/
```

不包含：第三方依赖、本机 profile、外部大网格、二进制、日志、结果和未提交工作树修改。

正式继续开发应从 `repository.bundle` clone；`source/` 更适合快速查看/编译：

```bash
git clone repository.bundle MPMC_SCW_case
git status --short --branch
```

## 11. 稳定点检查

内部稳定点前至少确认：

```text
[ ] git status 符合预期
[ ] make audit PASS
[ ] 没有机器私有路径/生成结果进入 Git
[ ] 受影响测试 PASS
[ ] PETSc-facing 改动已在目标环境验证
[ ] representative case 已记录
[ ] docs/CHANGELOG.md 只在形成重要能力/架构里程碑时更新
[ ] handoff package 能通过 manifest/SHA 验证
```
