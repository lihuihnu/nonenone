# AdaptiveTimeStepper/Core 子模块说明书

纯 C++ 时间步控制核心。外部 API：`AdaptiveTimeStepConfig`、`AdaptiveTimeStepPolicy`、`AdaptiveTimeStepper`、solve/event/statistics types。内部 API：retry、dt clipping/growth、fixed-target advance。Core 只依赖 backend contract，不知道 PETSc/Natural/井。任何新策略必须保持 begin-attempt / accept / reject 的事务边界。验证见 `adaptive_timestep_test.cpp`。
