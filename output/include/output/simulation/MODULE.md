# Output/Simulation 子模块说明书

负责模拟级派生诊断。外部 API 为 `ModelInventoryOutput`、`ComponentMassBalanceOutput`、`ReservoirDiagnosticsOutput` 及各自 Options。它们消费 runtime 提供的 global inventory/diagnostics，不参与求解。内部 API 是列定义、sample scheduling 和 CSV file state。新增诊断优先复用 runtime 已有 global query，避免再次遍历全网格/重算物性。验证见 `output_feature_test.cpp`。
