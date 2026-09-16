# AdaptiveTimeStepper/PETSc 子模块说明书

负责把一次 PETSc SNES 求解转换为 `NonlinearSolveResult`。外部 API 为 `PetscSnesDriver`。内部 API 是 SNES reason、nonlinear/linear iteration 统计和 Newton record 采集。它不决定是否 cut/grow dt；决策仍属于 core policy。改变 SNES 统计解释时必须同步 solver diagnostics。最终验证需要 PETSc integration。
