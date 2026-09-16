# MPMC_SCW Test

`test/` 是项目唯一正式测试入口。业务模块不再各自维护独立测试 Makefile。

## 1. 结构

```text
test/
├── Makefile
├── run.sh
├── src/
│   ├── unit/          PETSc-free 确定性单元/回归
│   ├── compile/       模板与配置实例化门禁
│   ├── integration/   PETSc/MPI/Grid 联调
│   ├── benchmark/     固定 workload 性能基准
│   └── support/       测试专用 ABI/fixture
└── include/test/
```

## 2. 测试层级

### Unit

覆盖 AD、common、indices、EOS/Flash、物性、active-set、质量守恒、井、输出以及代表性 case preflight。

```bash
make unit
```

或：

```bash
cd test
make run-unit CXX=g++
make run-unit CXX=clang++
```

### Compile matrix

只验证重要模板/模型组合能完整实例化，不把“能编译”当成数值正确。

### Integration

验证 PETSc/MPI/METIS/Grid、并行通信、装配和真实数据路径：

```bash
make integration NP=1 MESH_DIR=/path/to/DQ_data
make integration NP=2 MESH_DIR=/path/to/DQ_data
make integration NP=4 MESH_DIR=/path/to/DQ_data
```

### Benchmark

用于观察固定 workload 的性能变化，不作为单独物理正确性证据。

## 3. 完整命令

根目录：

```bash
make audit
make unit
make cases
make integration NP=2 MESH_DIR=/path/to/DQ_data
make full NP=1 MESH_DIR=/path/to/DQ_data
```

测试目录：

```bash
cd test
make clean
make run-unit CXX=g++
make full -j
make run-full NP=4 MESH_DIR=/path/to/DQ_data
make benchmark-performance CXX=g++
```

## 4. 新增测试规则

- 纯数学/物理逻辑优先写 unit；
- 新配置组合至少有 compile gate；
- MPI/ghost/partition/PETSc 行为必须写 integration；
- EOS/Flash 回归同时考虑公开数据或 tools cross-check；
- bug 修复优先先写能复现旧问题的最小回归；
- 测试不复制 production 公式作为“参考实现”，除非明确是独立算法交叉验证。

## 5. 验证结论怎么记录

不要只写“PASS”。至少记录：

```text
commit
compiler
PETSc/MPI
rank count
dataset/case
command
numerical tolerance / key metric
result
```

文档中的历史 `33/33` 等数字只是某一提交/环境的记录；当前状态必须重新运行测试确认。

## 6. 并行日志

Slurm/PMI/MPI launcher 可能输出与模型无关的环境信息。过滤日志时只允许去除已知 launcher 噪声，不能过滤 SNES divergence、MPI error、mass-balance error 或程序返回码。

详细测试职责见 [MODULE.md](MODULE.md)，机器配置见 [../docs/HPC_RUN.md](../docs/HPC_RUN.md)。
