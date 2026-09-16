# DQcase

> **定位**：真实/外部网格生产路径基线，重点覆盖 GRDECL/CpGrid、METIS、input-cell 编号、多井控制和 Legacy 三相多组分物理。

## 模型摘要

- 8 个 EOS 组分：`CO2, N2+C1, C2-nC4, iC5+nC5+, C7-C17, C18-C22, C23-C27, C28-C80`，另有独立水相。
- EOS：PR；Legacy aqueous-CO2 dissolution 与 Land trapping 启用。
- 初始状态：`22.91437 MPa`、`371.7611 K`，`Sw=0.55`、`So=0.45`、`Sg=0`。
- 网格从 `-mesh_dir` 或 `-grdecl` 获取；生产路径只由 rank 0 读取/构造/分区，非 root 接收 owned+one-ring ghost snapshot；`cell_id_map` 等一次性全局元数据写出后 root 也压缩为 local snapshot。direct GRDECL 在 Mesh 构造后只保留轻量 PORO/PERM payload，避免 expanded corner geometry 长期驻留。机器私有路径只放 `config/hpc.local.mk`。
- 完井 cell 使用**原始输入的 0 基 cell id**，运行时再映射到 CpGrid current id。

## 权威来源

`case_config.hpp` 保存流体/初值/时间；`well_config.hpp` 保存当前 DQ 对齐井表、WI、rate/BHP limits 与 schedule。README 不复制完整井表。

## 运行

配置好 DQ 数据路径后，从仓库根目录：

```bash
make case CASE=DQcase -j
make run CASE=DQcase NP=4
```

也可直接给 executable 传 `-mesh_dir` 或 `-grdecl`。环境与 CpGrid 数据配置见 `docs/HPC_RUN.md`。

## 验收

先运行 `legacy_dq_case_config_test.cpp`，再用 `make distributed-mesh NP=4` 验证自包含 root-only ingest/scatter；正式验证重点检查 root-only GRDECL/rock 读取、input-order 输出、多井控制切换、Land/dissolution 历史事务和 MPI 一致性。
