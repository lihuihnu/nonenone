# CpGrid 子模块说明书

## 1. 职责

`cpgrid` 实现角点/非结构 reservoir 网格：GRDECL 读取、网格拓扑、几何、rock field、METIS partition、owned/ghost DOF 及 PETSc Vec/Mat 支撑。

## 2. 外部 API

`cpgrid.hpp`（`CpGridCore` 生产入口，`CpGrid<Tag>` 兼容层）、`distributed_mesh_loader.hpp`、`distributed_rock_loader.hpp`、`cartesian_cell_directory.hpp`、`grdecl.hpp`、`rock_loader.hpp`、`dof_map.hpp`、`metis_partitioner.hpp`、`petsc_io.hpp`。

### `CpGridCore`

主要能力：

- 装载/建立 `Mesh` 并完成 `setup()`。
- `registerLayout(dofPerCell)` 注册额外 per-cell layout。
- 创建 global/local PETSc vector 与 matrix 所需布局。
- 访问 owned/local snapshot cells、faces、cell geometry、rock property。
- `rockLocalView()`：RAII 同步 owned+ghost 孔隙度/渗透率。
- 输入编号与 current/partition 编号转换。

### GRDECL API

`grdecl.hpp`：

- `GrdeclUnitSystem`
- `GrdeclLoadOptions`
- `GrdeclGridData` / `GrdeclCell` / `GrdeclPoint`
- GRDECL loader 将 `SPECGRID/DIMENS`、`COORD`、`ZCORN`、`ACTNUM`、`PORO/PERM*` 转为 canonical 数据。

`GrdeclLoadOptions::keepLargestConnectedComponent` 默认为 `false`，因此通用读取不会静默删除
`ACTNUM` 单元。显式启用后，loader 按 Cartesian 六邻域保留最大的面连通分量，并在
`GrdeclGridData::disconnectedCellCountRemoved` 报告删除数量；同规模分量使用最小 Cartesian
编号确定结果。DQ direct-GRDECL 路径启用该策略，以复现 MRST
`buildGridFromGRDECL(..., 'largest')` 的网格选择。

`CpGridRockLoadOptions` 的零值替换是显式输入清洗策略。DQ direct-GRDECL 路径对孔隙度及三个
方向渗透率使用正值样本均值替换零值，与参考 `DQ_data` 的生成过程一致；已经归一化的 CSV
路径不会再次修改属性。

### DOF/partition API

- `DofLayout`、`DofMap`
- `LayoutRegistry`
- `MetisPartitioner`
- `Sparsity`

## 3. 内部 API

`Mesh`、`Polyhedron`、`Face`、`Node`、`Point` 及 `grid/src/*.cpp` 中的拓扑构造属于内部实现。Natural 层应通过 backend 访问，不应依赖 Mesh 内部容器布局。

## 4. 整体逻辑

```text
root-only CSV / GRDECL
  -> tokenizer / parser / unit and MAPAXES normalization
  -> ACTNUM filtering / corner-point topology builder
  -> Mesh cells/faces/nodes
  -> root topology / METIS partition
  -> owned + ghost MPI snapshot
  -> all-rank local topology
  -> DofMap
  -> rock vectors
  -> cached pure geometry / TPFA geometric factor
  -> NaturalCpGridBackend
```

GRDECL corner coordinates 由 `COORD` pillar 与 `ZCORN` 深度插值得到；`ACTNUM=0` cell 不进入 active grid。

## 5. 扩展规则

- 新输入 keyword 应首先进入 `grdecl.cpp` parser/normalizer。
- `CpGridCore` 自身只注册 rock 需要的 3/1 DOF；Natural 主变量与 phase-state 由 backend/runtime 显式注册。
- 新 rock field 不应改变已有 porosity/permeability API；需要时新增独立 layout。
- 不要把 Eclipse 完整语义（例如 NNC、LGR）假装成当前 conforming logical-neighbor 拓扑；若实现必须同步更新说明和测试。

## 6. 不可破坏约束

- `createInputOrderedCopy()` 对外保持输入 active-cell 顺序。
- owned/ghost 数组必须通过 RAII 或严格 restore。
- topology/geometry 在 `setup()` 后视为静态，Natural runtime 可安全缓存其派生量。
- DQ 在一次性全局编号表写出后，root 也压缩为 owned+ghost；direct GRDECL 只保留轻量 PORO/PERM payload，expanded corner geometry 必须尽早释放。
- 全局 Cartesian 目录使用连续只读数组和二分查找；单调 input/Cartesian 顺序时不保留额外 permutation。
- compact snapshot 的 current/input/node 本地定位使用 owned 连续区间及排序 ghost/node；DofMap ghost 查询同样使用排序数组。
- rank-major current-id 前缀由 Mesh 在分区后缓存并随 distributed directory 广播，DofMap 不为每种 DOF layout 重复收集。
- CpGrid setup 使用融合 cell/face 几何内核，一次 traversal 同时得到 volume/centroid 或 area/centroid/normal；历史独立几何 API 保留并由 bit-level 单测锁定等价性。
- `DofLayout` 注册阶段只构造实际需要的 global-to-local scatter；Jacobian `Sparsity` 与 local-to-global scatter 均按首次使用延迟创建。

## 7. 验证

`test/src/unit/grdecl_test.cpp`、`grdecl_rock_payload_test.cpp`；PETSc/METIS 行为由 `make distributed-mesh NP=4` 和 `test/src/integration` 覆盖。
