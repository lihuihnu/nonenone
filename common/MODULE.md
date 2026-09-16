# Common 模块说明书

## 1. 职责

`common` 保存不带 reservoir 物理语义的基础设施：数学工具、插值、区间搜索、单位、控制台输出、小容量容器、类型 traits，以及少量统一 PETSc 文件 I/O。

## 2. 外部 API

主要头文件：

- `common/math.hpp`：`MathToolbox`、通用数学辅助函数。
- `common/piecewise_linear_interpolation.hpp`：`PiecewiseLinearInterpolation`。
- `common/bilinear_interpolation.hpp`：`BilinearInterpolation`。
- `common/interval_search.hpp`：区间定位与 `BoundaryPolicy`。
- `common/units.hpp`：工程单位换算常量/函数。
- `common/small_vector.hpp`：`SmallVector`。
- `common/console.hpp`：`ConsoleSection` 与统一控制台排版。
- `common/petsc_io.hpp`：PETSc 向量读写基础操作。

## 3. 内部 API

- 插值器内部的区间查找与边界处理。
- `math.hpp` 中的返回类型推导辅助模板。
- `petsc_io.cpp` 中对 PETSc Viewer/Vec 的具体调用。

业务模块不应在 `common` 中加入“气相、井、flash、组分”等模型特定概念。

## 4. 整体逻辑

`common` 是横向基础层：AD、Grid、Natural、Output、Case 可以使用它，但它本身不反向依赖这些模块。

## 5. 扩展规则

只有满足“至少两个上层模块都可复用、且不包含领域语义”的工具才适合进入 `common`。单个模块私有 helper 应留在原模块。

## 6. 不可破坏约束

- 单位 API 必须明确输入/输出单位。
- 插值器不得暗中改变边界策略。
- PETSc I/O helper 必须成对恢复 PETSc 数组/Viewer 资源。

## 7. 验证

`test/src/unit/common_test.cpp` 和 `ad_common_integration_test.cpp`。
