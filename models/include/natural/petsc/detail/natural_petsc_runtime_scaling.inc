/**
 * @file natural_petsc_runtime_scaling.inc
 * @brief Natural/PETSc 方程缩放和变量缩放实现。
 *
 * 本文件是类体内实现片段，只能由对应公共头文件 include，不应独立包含。
 */

    [[nodiscard]] double variableScaleFactor_(int variable) const noexcept
    {
        if (!options_.scaling.enabled)
            return 1.0;
        if (variable == Indices::Primary::pressure ||
            (Indices::hasWellUnknown && variable == Indices::Primary::wellPressure))
        {
            return options_.scaling.pressureScale;
        }
        if (variable == Indices::Primary::liquidSaturation ||
            variable == Indices::Primary::vaporSaturation ||
            (Indices::hasWater && variable == Indices::Primary::waterSaturation))
        {
            return options_.scaling.saturationScale;
        }
        return options_.scaling.compositionScale;
    }

    [[nodiscard]] double equationScaleFactor_(
        int equation,
        PetscInt localBlock) const noexcept
    {
        if (!options_.scaling.enabled)
            return 1.0;

        if (equation >= Indices::Equation::massConservationBegin &&
            equation < Indices::Equation::afterMassConservation)
        {
            return 1.0 / options_.scaling.massResidualScale;
        }
        if constexpr (Indices::hasIndependentWaterConservation)
        {
            if (equation == Indices::Equation::waterConservation)
                return 1.0 / options_.scaling.massResidualScale;
        }

        if (equation >= Indices::Equation::fugacityBegin &&
            equation < Indices::Equation::fugacityBegin + Indices::numComponents)
        {
            return 1.0 / options_.scaling.fugacityResidualScale;
        }
        if constexpr (Indices::fullyCompositionalThreePhase)
        {
            if (equation >= Indices::Equation::waterFugacityBegin &&
                equation < Indices::Equation::waterFugacityBegin + Indices::numComponents)
            {
                return 1.0 / options_.scaling.fugacityResidualScale;
            }
        }
        if constexpr (Indices::hasAqueousCO2Dissolution)
        {
            if (equation == Indices::Equation::aqueousCO2Equilibrium)
                return 1.0 / options_.scaling.fugacityResidualScale;
        }
        if (equation == Indices::Equation::volumeClosure)
            return 1.0 / options_.scaling.closureResidualScale;

        if constexpr (Indices::hasWellUnknown)
        {
            if (equation == Indices::Equation::wellControl)
            {
                const std::size_t block = static_cast<std::size_t>(localBlock);
                if (block < representativeWellByLocalBlock_.size())
                {
                    const std::size_t wellIndex =
                        representativeWellByLocalBlock_[block];
                    if (wellIndex != std::numeric_limits<std::size_t>::max())
                    {
                        const Well &well = wells_[wellIndex];
                        if (well.isActive(currentTime_) && well.isRateControlled())
                        {
                            return 1.0 / std::max(
                                std::abs(well.signedTarget()),
                                options_.scaling.rateWellResidualFloor);
                        }
                    }
                }
                return 1.0 / options_.scaling.pressureScale;
            }
        }
        return 1.0;
    }

    // G8I 直接在 self/face/well 的 Jacobian pack 路径应用行、列缩放，
    // 数学上等价于原左右对角缩放，同时避免两个辅助 PETSc Vec 和一次
    // 完整稀疏矩阵 MatDiagonalScale 遍历。
    /** @brief 任意 rank 都可计算某口井控制行的缩放，不依赖代表单元本地化。 */
    [[nodiscard]] double wellControlEquationScaleFactor_(
        std::size_t wellIndex) const noexcept
    {
        if (!options_.scaling.enabled)
            return 1.0;
        if (wellIndex >= wells_.size())
            return 1.0 / options_.scaling.pressureScale;

        const Well &well = wells_[wellIndex];
        if (well.isActive(currentTime_) && well.isRateControlled())
        {
            return 1.0 / std::max(
                std::abs(well.signedTarget()),
                options_.scaling.rateWellResidualFloor);
        }
        return 1.0 / options_.scaling.pressureScale;
    }
