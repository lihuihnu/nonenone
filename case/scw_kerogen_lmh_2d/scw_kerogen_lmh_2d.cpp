/**
 * @file scw_kerogen_lmh_2d.cpp
 * @brief 水–轻中重烃二维三状态方程对比算例入口。
 */
#include <case/petsc_custom_hooks.hpp>
#include <case/structured_three_eos_compare_runner.hpp>

#include "case_config.hpp"
#include "well_config.hpp"

namespace Case
{
struct Definition final
{
    using Config = ScwKerogenLmh2D::Config;
    using PrFactoryConfig = ScwKerogenLmh2D::PrFactoryConfig;
    using SwFactoryConfig = ScwKerogenLmh2D::SwFactoryConfig;
    using CpaFactoryConfig = ScwKerogenLmh2D::CpaFactoryConfig;
    [[nodiscard]] static constexpr const auto &wells() noexcept
    {
        return WellConfig::wells;
    }
};
using Runner = MPMC::cases::StructuredThreeEosCompareRunner<Definition>;
using Runtime = typename Runner::Runtime;
} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    return Case::Runner::runPetscMain(argc, argv);
}
