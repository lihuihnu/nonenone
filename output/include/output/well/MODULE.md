# Output/Well 子模块说明书

负责井历史和井控切换输出。外部 API 为 `DetailedWellOutput`、`WellControlSwitchOutput`、`WellHistory` 及 Options。内部 API 是列布局、单位换算和文件状态。数据来自 `WellState`/control events；不得自行重新评价井方程。注入/生产符号与工程 magnitude 的转换必须在列名和文档中明确。验证见 `output_well_test.cpp`。
