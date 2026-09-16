# StructuredGrid 子模块说明书

## 1. 职责

`structuredgrid` 提供规则三维笛卡尔网格和 PETSc DMDA 资源，负责几何、activity、rock property、六邻接面及结构化 global/local 索引。

## 2. 外部 API

- `StructuredGridCore`：生产 façade，负责跨组件编排、active/face/TPFA 与资源收尾；`StructuredGrid<Tag>` 仅作兼容层。
- `GridExtent`：按物理范围指定网格。
- `CellSize`：按均匀单元尺寸指定网格。
- `StructuredFace` / `StructuredFaceOrder`：六面邻接信息。
- `StructuredGridIndex`：`i,j,k` 与 global cell id 转换。
- `DMDAContainer`：DM 生命周期和 vector/matrix 相关资源。
- `DMDAElem`、`DMDARegion`、`DMDADimensions`：本地/全局 DMDA 区域描述。
- `structuredgrid_components.hpp`：轴向几何状态、DMDA/PETSc layout、active 与岩石 Vec 的内部职责组件。
- `structured_cartesian_geometry.hpp`：PETSc-free 的体积、面半程距离/面积和轴中心计算。
- `structuredgrid_views.hpp`：岩石属性读写 RAII 视图与装配期数组访问守卫。

## 调用关系

`dimensions/extent -> Grid canonical DMDA -> active ghost state + derived geometry -> rock views -> NaturalStructuredGridBackend -> model layout registration`。

`StructuredGridCore` 对上层暴露 setup、cell/face 查询、geometry、porosity/permeability、DM/Vec/Mat 支撑以及 activity 信息。

Natural 装配应通过 `AssemblyAccess` 管理 ghosted 岩石数组生命周期；体积、面面积和半程距离由轴宽按需计算，重力只保留 O(nz) 垂向中心缓存。case 初始化岩石优先使用 `rockWriteView()`，报告/诊断优先使用 `rockReadView()`。

## 3. 内部 API

几何 Vec 建立、DMDA array access、local cache、resource release、activity 初始化等私有函数属于实现细节。

- `StructuredGridCore::setup()` 只注册 Grid 自己的 1/3 DOF 布局；6/18 DOF 仅由旧 `StructuredGrid<Tag>` 兼容 façade 为历史 geometry Vec 注册。Natural 主变量、AD 辅助布局和 phase-state 由 backend/runtime 显式注册。
- 项目内部新代码不得直接依赖 `permeabilityVector`、`porosityVector`、`dx/dy/dz`、`Lx_/Ly_/Lz_` 或手工 `mapLocalArrays()/unmapLocalArrays()`；旧 geometry Vec/raw arrays 只保留在 `StructuredGrid<Tag>` façade，不再进入生产 Core。
- 普通求解向量仍可通过通用 `vecGetArray()/vecRestoreArray()` 访问；不要把岩石属性重新退回裸 Vec 管理。
- 六面顺序统一由 `StructuredFaceOrder` 表达；Cartesian global id 规则保持稳定。

## 4. 整体逻辑

```text
nx,ny,nz + extent/cell-size
  -> DMDA
  -> cell centers / volumes / face normals / distances
  -> activity
  -> porosity/permeability vectors
  -> local ghost snapshots
  -> NaturalStructuredGridBackend
```

## 5. 扩展规则

结构化网格新增 property 时优先增加 per-cell Vec，而不是在 Natural 中按 `(i,j,k)` 特判。六面顺序必须继续使用 `StructuredFaceOrder` 统一定义。

## 6. 不可破坏约束

- DMDA array 获取和恢复必须配对。
- global index 规则不能随意改变，否则输出和井 cell id 会失配。
- 几何初始化后不得在 Newton 中修改。

## 7. 验证

StructuredGrid 的核心行为由 Natural unit/integration case 间接覆盖；完整 PETSc 行为需运行 `structured_grid_test`、`natural_structured_test` 及真实 PETSc/MPI case 验证。
