# Models 模块说明书

`models` 是工程领域模型层，目前包含 `natural` 和 `well` 两个相互解耦的主模块。`well` 定义井的纯数据/控制模型，`natural` 通过 `physics/well_source.hpp` 和 PETSc runtime 将井耦合进储层方程。外部代码应分别使用 `natural/natural.hpp`、`natural/petsc/natural_petsc.hpp` 和 `well/well.hpp`；不要新建跨越两者职责的第三套井/流体实现。详细 API 见 `include/natural/MODULE.md` 与 `include/well/MODULE.md`。
