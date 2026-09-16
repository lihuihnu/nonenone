# Output/VTK 子模块说明书

负责可选 CpGrid VTK 导出。外部 API 为 `CpGridVtkWriter` 与 `VtkCellOrdering`。内部 API 是 VTK connectivity、cell ordering 和 field serialization。VTK 是可选表现层，不得成为求解或标准 CSV 的依赖；cell ordering 必须显式选择并与 input/current id 语义一致。最终验证需要带 VTK 的 integration build。
