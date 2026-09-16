# AD 模块说明书

## 1. 职责

`ad` 提供 Natural 全隐式离散使用的前向自动微分标量。它只负责“数值 + 一阶导数”的代数传播，不知道组分、网格、PETSc 或井。

## 2. 外部 API

推荐直接包含：

```cpp
#include <ad/Evaluation.hpp>
#include <ad/DynamicEvaluation.hpp>
#include <ad/Math.hpp>
```

主要类型：

- `DenseAd::Evaluation<ValueT, numVars, staticSize>`：固定导数维度 AD；Natural 主路径使用该类型。
- `DenseAd::Evaluation<ValueT, DynamicSize, staticSize>`：运行期导数维度版本。
- `DenseAd::MathToolbox<...>`：为 AD 标量提供数学适配。
- `DenseAd::is_evaluation<T>`：类型判定。

常用操作包括 `createVariable(value, derivativeIndex)`、`value()`、`derivative(i)`、`setDerivative()`、`clearDerivatives()` 以及普通算术/数学函数重载。

## 3. 内部 API

`Math.hpp` 内的数学函数转发、返回类型推导和 `MathToolbox` 特化属于实现细节。业务模块不应依赖内部辅助模板名称，只依赖标准数学表达式可用于 `Evaluation`。

## 4. 整体逻辑

```text
double primary value
  -> Evaluation::createVariable
  -> EOS / properties / flux / accumulation 使用普通表达式
  -> 运算符自动传播导数
  -> residual.value + residual.derivative(j)
  -> PETSc residual / Jacobian
```

固定维度模式使每单元 Jacobian block 的导数布局在编译期确定；动态维度主要服务通用测试和非固定尺寸场景。

## 5. 扩展规则

- 新数学函数应在 `Math.hpp` 中实现 AD 对应重载，而不是在物理模块中手工拆 value/derivative。
- 不要在物理路径中通过 `double(...)` 或 `scalarValue()` 提前丢失导数，除非该量按数学定义必须被冻结。
- 改变导数维度前必须同步检查 `Indices::ValueType` 和 Jacobian 装配。

## 6. 不可破坏约束

- `Evaluation<double,N>` 的算术结果必须与 double 路径数值一致。
- 导数索引必须与 `Indices::Primary` 一一对应。
- 不允许 AD 模块依赖 PETSc。

## 7. 验证

`test/src/unit/ad_test.cpp`、`ad_dimension_test.cpp`、`ad_common_integration_test.cpp`。
