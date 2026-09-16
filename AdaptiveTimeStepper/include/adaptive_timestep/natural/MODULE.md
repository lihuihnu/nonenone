# AdaptiveTimeStepper/Natural 子模块说明书

把 Natural runtime 映射到通用 stepper backend contract。外部 API 为 `NaturalAdaptiveBackend<Runtime,...>`、`NoAcceptedStepHook` 与 `NoFailedSolveHook`。它转发 current time、dt、solve、accept/reject、well-control 更新，并在 accepted step 触发可选 hook。内部 API 是 runtime/well-cycle 协调。不得在此复制 adaptive policy 或 Natural 物理。验证由 adaptive/well tests 和正式 case 覆盖。

失败求解 hook 必须在 `rejectAttempt()` 之前执行，只允许读取失败试探状态做诊断，不得提交历史量或修改 adaptive policy。
