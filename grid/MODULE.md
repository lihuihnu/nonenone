# Grid 模块说明书

## 1. 职责

`grid` 只负责空间离散数据：拓扑、几何、岩石属性、单元/面编号、DOF 映射、并行 partition 和 PETSc 网格资源。它不负责 EOS、相平衡、井控制或时间推进。

## 2. 外部 API

工程提供两套 backend：

- `CpGridCore`：角点/非结构生产路径；`CpGrid<Tag>` 仅保留旧源码兼容，见 `grid/include/cpgrid/MODULE.md`。
- `StructuredGridCore`：规则笛卡尔/DMDA 生产路径；`StructuredGrid<Tag>` 仅保留旧源码兼容，见 `grid/include/structuredgrid/MODULE.md`。

生产 CpGrid 路径为：`root-only MRST/GRDECL -> root topology/METIS/current-id -> owned+ghost MPI snapshot -> DOF -> Natural backend`。

旧 `Mesh::prepareForUse()` replicated 路径继续保留用于兼容与 A/B 验证。

上层 Natural/PETSc 不直接依赖两者具体实现，而通过 `NaturalCpGridBackend` / `NaturalStructuredGridBackend` 适配。

## 3. 内部 API

新输入格式先转换为 `grid/src` 内部的 `CanonicalMeshData`，再物化为 Mesh；不要在 Mesh topology 或 Natural 中增加格式分支。

`grid/src` 中的 Mesh、Polyhedron、Face、DOF、METIS、PETSc resource 实现属于网格内部。除 backend 外，上层不应直接拼装这些内部对象。

## 4. 整体逻辑

```text
输入网格
  -> topology + geometry
  -> active cells
  -> porosity/permeability
  -> partition / ghost
  -> DOF layout
  -> face connection + transmissibility geometry
  -> Natural backend 消费
```

## 5. 扩展规则

新增网格格式应先转换为已有 Grid 的 canonical cell/face/rock 表示，避免在 Natural 中加入格式分支。新增网格 backend 必须提供与当前 backend 等价的 cell iteration、local/global DOF、geometry、rock、vector/matrix 构造能力。

## 6. 不可破坏约束

- 外部结果编号使用输入/canonical cell id，而不是 MPI current id。
- 网格层不理解 phase/component。
- 几何和岩石属性单位必须进入 Natural 前统一。

## 7. 验证

GRDECL 由 `grdecl_test.cpp` 覆盖；PETSc/MPI 网格行为由 integration tests 覆盖。GridCore 不理解 phase/component/AD，也不读取 Natural Tag；模型布局由 backend 显式注册。`make distributed-mesh NP=4` 不依赖外部 DQ 数据，可直接验证 root-only ingest/scatter。
