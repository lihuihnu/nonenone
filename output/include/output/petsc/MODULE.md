# Output/PETSc 子模块说明书

负责从 PETSc/CpGrid 分布式向量中安全收集输出数据。外部 API 为 `CpGridOutputAccess`、`CpGridSaver`、`OutputPetscOwnedReadView`。内部 API 包括 owned array view、input/current id remap 和 MPI collection。主逻辑是“distributed current-order Vec -> owned read -> canonical/input order -> writer”。不得在此修改 Vec 或计算 EOS；所有数组访问必须 RAII restore。最终验证需要 PETSc/MPI integration。
