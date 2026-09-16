# AdaptiveTimeStepper/Well 子模块说明书

负责“同一 timestep 内井控切换并重新求解”的事务。外部 API 为 `WellControlCycle`、`NoWellControlCycle`、`WellControlSnapshot` 和 switch callback。内部 API 是 control snapshot/restore、limit evaluation 和 switch notification。井控尚未 settled 时不能接受 timestep；rejected attempt 必须恢复原 control。验证见 `adaptive_well_test.cpp`。
