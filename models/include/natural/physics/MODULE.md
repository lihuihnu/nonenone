# Natural Physics 子模块说明书

## 1. 职责

`physics` 保存进入守恒方程或 constitutive state 的独立物理过程。

## 2. 外部 API

- `accumulation.hpp`：`AccumulationResult` 与当前/历史质量库存计算。
- `face_flux.hpp`：`FaceMassFlux`、`UpwindSide` 和 TPFA + upwind 组分通量。
- `well_source.hpp`：`PerforationWellResult` 与 perforation source。
- `aqueous_co2.hpp`：Legacy `AqueousCO2Model`。
- `competitive_adsorption.hpp`：`CompetitiveLangmuirAdsorption<N>`。
- `land_trapping.hpp`：`LandTrappingModel`。

## 3. 内部 API

上游选择、phase contribution 汇总、Land free/trapped saturation 计算等 helper 属于内部实现。外部不应复制这些公式。

## 4. 整体逻辑

```text
CellState + CellProperties
  -> accumulation
Neighbor properties + face geometry
  -> phase potential -> upwind -> component flux
Well + perforation + local properties
  -> well source
Optional model flags
  -> Land / dissolution / adsorption contribution
```

## 5. 扩展规则

新增过程若改变守恒库存/通量/源汇，应在本层建立独立可测试函数，再由 `CellPropertyEvaluator` 或 assembly 调用。Fully-compositional EOS 已处理 H2O/CO2 相间分配，不得再叠加 legacy aqueous-CO2 equilibrium。

## 6. 不可破坏约束

- 井/通量符号约定保持“注入为正、生产为负”。
- TPFA 上游判断与原公式保持一致。
- Land 历史量只在 accepted timestep 更新。

## 7. 验证

Natural extended/conservation、well_natural、Land/溶解相关 case tests。
