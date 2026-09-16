# Test 模块说明书

## 1. 职责

`test` 是工程唯一正式测试入口。测试分为纯 C++ unit、模板/配置 compile matrix、PETSc/MPI integration 和 performance benchmark。

## 2. 外部使用 API

常用命令：

```bash
cd test
make clean
make run-unit CXX=g++
make run-unit CXX=clang++
make full -j
make run-full NP=4 MESH_DIR=/path/to/mesh
make benchmark-performance CXX=g++
```

目录：

- `src/unit/`：确定性纯 C++ 回归。
- `src/compile/`：只验证模板/配置能完整实例化。
- `src/integration/`：真实 PETSc/MPI/Grid 运行。
- `src/benchmark/`：固定 workload 局部性能回归。
- `src/support/`：测试共享数据/构造器。

## 3. 内部 API

测试 helper 只能服务测试，不得反向被 production include。`test/include/test` 中的工具不是生产 API。

## 4. 测试策略

```text
发现 bug
  -> 最小 deterministic reproduction
  -> unit regression
  -> production fix
  -> related unit
  -> all unit clean build
  -> compile matrix
  -> integration / formal case
```

物理模型还需要 literature/reference validation，而不仅是“代码不崩”。

## 5. 新增测试规则

- 新算法优先 unit。
- 新 compile-time feature combination 增加 compile test。
- 必须依赖 PETSc/MPI 的行为放 integration。
- performance test 必须使用固定 workload + checksum，不能以改变结果换速度。

## 6. 不可破坏约束

- 不为了测试方便污染 production API。
- header-only 重大修改必须 clean rebuild，防止 stale binary 假通过。
- `make full` 只代表完整编译；不能冒充 `run-full` 已运行。

## 7. 当前基线

当前具有 33 个纯 C++ unit executables；除基础相图与 bubble/dew 外，还覆盖 O/G/W phase role、SW fail-only recovery、K-value、critical-estimate 质量字段、经验参数、EOS regression，以及五组分独立优选参数算例的配置/闪蒸预检。内部稳定点仍要求 clean build 全部通过。
