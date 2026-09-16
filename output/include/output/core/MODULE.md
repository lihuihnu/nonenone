# Output/Core 子模块说明书

负责文件和 CSV 基础设施，不理解 Natural/PETSc。外部 API 为 `OutputWriter`、`OutputTextFile`、`OutputCellSelection`、`OutputCellIdSpace`、`MassSeries` 及 `format.hpp` 的格式化函数。内部 API 是 header 状态、stream 生命周期和 CSV row/header helper。调用链为“typed values -> format/header -> stream”。新增格式化能力应保持无模型语义；不得在本层做 MPI reduction 或重算物性。验证见 `output_test.cpp`。
