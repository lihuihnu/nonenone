# 本地 Linux / WSL 与 HPC 编译运行

本项目只维护一套 Makefile 工作流。机器差异通过 `config/hpc.mk` 与私有 `config/hpc.local.mk` 处理，不把个人绝对路径写进源码。

## 1. 配置层

默认配置：

```text
config/hpc.mk
```

本机覆盖：

```text
config/hpc.local.mk
```

首次配置：

```bash
cp config/local.example.mk config/hpc.local.mk
```

常见字段：

```makefile
MPMC_PLATFORM ?= local
PETSC_DIR      ?= /path/to/petsc-3.22.2
PETSC_ARCH     ?= arch-linux-c-opt
FMT_DIR        ?= /path/to/fmt
METIS_DIR      ?= /path/to/metis
MPI_LAUNCHER   ?= mpiexec
MPI_LAUNCHER_ARGS ?=
MODULES_ENABLED ?= 0
DQ_MESH_DIR    ?= /path/to/DQ_data
```

检查最终解析值：

```bash
make print-config
make doctor-local
# 或
make doctor-hpc
```

`config/hpc.local.mk` 不提交 Git。

## 2. 依赖

基础：

- C++17 compiler；
- GNU Make；
- Python 3；
- PETSc 3.22.x；
- 与 PETSc ABI 一致的 MPI；
- METIS；
- fmt headers。

VTK/X11 只在相关完整门禁需要时配置。

Windows 推荐 WSL2/Ubuntu；不维护原生 MSVC + PETSc 路径。

## 3. 不依赖 PETSc 的检查

```bash
make audit
make unit
```

需要分别验证 GCC / Clang 时：

```bash
make unit-gcc
make unit-clang
```

## 4. 编译正式算例

全部：

```bash
make cases -j
```

单个：

```bash
make case CASE=three_eos_3d_compare -j
```

`case/Makefile` 是正式算例唯一构建入口；不要在算例目录再维护第二套依赖配置。

## 5. 本地运行

推荐从根目录：

```bash
make run CASE=three_eos_3d_compare NP=2 EOS=pr RESULT_DIR=./results/pr
```

也可先生成算例自己的 managed `run.sh`：

```bash
make prepare CASE=three_eos_3d_compare TASKS=2
```

运行参数、环境变量和结果目录以 [../case/README.md](../case/README.md) 为准。

## 6. HPC 提交

```bash
make prepare CASE=three_eos_3d_compare TASKS=8
make submit  CASE=three_eos_3d_compare TASKS=8
```

launcher 与 scheduler 由 profile 决定。不要把某个集群专用参数硬编码进通用 `run.sh`。

特别注意：

- `mpiexec`、`mpirun`、`yhrun` 的参数不能混用；
- PETSc 必须与运行时 MPI ABI 匹配；
- `--mpi=pmi2` 等调度器参数不能传给普通 MPICH/OpenMPI `mpiexec`；
- module 系统只在 `MODULES_ENABLED=1` 时使用。

## 7. PETSc/MPI 集成门禁

```bash
make integration NP=1 MESH_DIR=/path/to/DQ_data
make integration NP=2 MESH_DIR=/path/to/DQ_data
```

改动网格、通信、装配或 ghost 逻辑时至少补：

```bash
make integration NP=4 MESH_DIR=/path/to/DQ_data
```

完整门禁：

```bash
make doctor-full
make full NP=1 MESH_DIR=/path/to/DQ_data
```

本地通过不能替代目标超算重复验证。

## 8. 常见故障

### Make 仍指向旧机器绝对路径

先检查：

```bash
make print-config
```

确认是否误提交/误加载了旧 `config/hpc.local.mk`。机器路径只应存在于 profile。

### `module: command not found`

本地配置：

```makefile
MODULES_ENABLED ?= 0
```

### `set: Illegal option -o pipefail`

脚本必须由 Bash 执行，不要用 `sh run.sh`。

### `mpiexec: command not found`

在 `config/hpc.local.mk` 设置正确 `MPI_LAUNCHER`，并确认它与 PETSc 使用的 MPI 一致。

### `mpiexec ... unrecognized argument mpi`

通常是把 Slurm/yhrun 参数传给了普通 MPI launcher。清空或修正 `MPI_LAUNCHER_ARGS`。

### `libpetsc.so: undefined reference to updateState/updateSol`

这是历史定制 PETSc hook ABI 问题。当前工程保留标准 PETSc callback 与历史 C/C++ hook 兼容路径；先用：

```bash
make doctor-local
```

确认链接到的是预期 PETSc，而不是系统中另一套库。

### 多进程才失败

优先检查：

1. PETSc/MPI ABI；
2. communicator；
3. owned/ghost；
4. partition/DOF；
5. collective 是否所有 rank 同步进入。

不要先通过放松 SNES/KSP 容差掩盖并行错误。

## 9. 清理

```bash
make clean
```

测试/工具也可在各自模块使用自己的 `make clean`。运行结果、日志和本机配置不进入 Git。
