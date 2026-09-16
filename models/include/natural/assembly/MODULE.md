# Natural Assembly 子模块说明书

## 1. 职责

`assembly` 将已经计算好的当前/上一时层单元物性、面通量和井源汇组织为固定大小的局部非线性方程向量。

## 2. 外部 API

- `CellResidualBlock<Indices,Scalar>`：一个 cell 的 `numEquations` 个 residual。
- `assembleCellResidual(...)`：统一局部组装入口。
- `local_equations.hpp`：fugacity equilibrium、inactive-phase replacement、saturation closure 等局部代数方程辅助。

## 3. 内部 API

方程索引选择、inactive phase 方程替换和残差缩放属于 assembly 内部；调用者只应通过 `Indices::Equation` 理解最终布局。

## 4. 整体逻辑

```text
accumulation(current,previous)/dt
 + face component flux sum
 - well source
 -> conservation equations

current phase fugacities / phase presence
 -> equilibrium or inactive-phase equations

saturations
 -> closure
```

## 5. 扩展规则

新增守恒量必须先在 `Indices` 增加 equation/primary 对应关系，再在 assembly 中加入唯一残差。不要仅为输出增加方程。

## 6. 不可破坏约束

局部系统维度固定；相消失不允许改变 PETSc block size，只能替换 inactive equation 内容。

## 7. 验证

Natural core/extended、two-cell conservation、three-phase AD tests。
