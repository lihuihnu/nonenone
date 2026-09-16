# Natural Properties 子模块说明书

## 1. 职责

`properties` 提供从 EOS/组分/压力到密度、质量分数、黏度等基础物性的模型封装。

## 2. 外部 API

- `BlackOilPropertyModel`：Legacy black-oil/compositional oil-gas 辅助物性。
- `CompositionalPropertyModel`：根据 `CompositionalMixture`/EOS 计算 fully-compositional 所需物性。
- `Iapws2008IndustrialAqueousViscosity`：Water 角色的 IAPWS R12-08
  工业黏度；非水摩尔分数超过显式适用域时拒绝计算，不作隐藏混合物外推。
- `McBrideWright2015AqueousViscosity`：受组成适用域保护的 H2O-CO2
  水相黏度闭包。

`FluidSystem` 持有这些模型，上层通常不直接构造第二份 property model。

## 3. 内部 API

临界性质组合、LBC 黏度常量与 composition conversion helper 属于实现细节。

## 4. 整体逻辑

```text
pressure + T + composition + EOS result
  -> molar density
  -> mass fractions
  -> mass density
  -> viscosity
```

## 5. 扩展规则

新增黏度/密度 correlation 应保持单位契约，并通过 `FluidSystem` 统一暴露。纯组分常数优先预计算，不要在每 Newton/cell 重算。

## 6. 不可破坏约束

物性模型不得决定 phase presence；相态由 state/thermo 层负责。若 Water
公共角色会选择专用水相物性，EOS/flash 必须同步配置同一组成适用域，防止
超出闭包适用范围的烃富相被误标为 Water。

## 7. 验证

Natural/thermo unit tests 以及 performance benchmark 的 cell-property workload。
