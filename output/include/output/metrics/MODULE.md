# Output/Metrics 子模块说明书

负责与文件格式无关的数值诊断。外部 API 为 `ComponentTotals`、`ComponentMassBalanceSnapshot`、`ComponentMassBalanceLedger`。Ledger 接受初始库存、accepted 时间段的注入/产出率和当前库存，计算累计量与闭合误差。内部 API 是 finite-check、time integration 和 initialized-state guard。只有 accepted timestep 才允许推进 ledger；rejected attempt 不能重复积分。验证见 `component_mass_balance_test.cpp` 与 output tests。
