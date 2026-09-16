/**
 * @file five_component_eos_tuned_compare.cpp
 * @brief 五组分独立优选参数 PR/SW/CPA 流动算例的运行时 EOS 选择入口。
 */
#include <case/petsc_custom_hooks.hpp>
#include "case_config.hpp"
#include "well_config.hpp"

#include <case/structured_three_eos_compare_runner.hpp>

namespace Case
{

struct Definition final
{
    using Config = CaseConfig::Config;
    using PrFactoryConfig = CaseConfig::PrFactoryConfig;
    using SwFactoryConfig = CaseConfig::SwFactoryConfig;
    using CpaFactoryConfig = CaseConfig::CpaFactoryConfig;

    [[nodiscard]] static constexpr const auto &wells() noexcept
    {
        return WellConfig::wells;
    }
};

using Runner = MPMC::cases::StructuredThreeEosCompareRunner<Definition>;
using Config = typename Runner::Config;
using ModelConfig = typename Runner::ModelConfig;
using Indices = typename Runner::Indices;
using Grid = typename Runner::Grid;
using Runtime = typename Runner::Runtime;
using Well = typename Runner::Well;

int run()
{
    return Runner::run();
}

} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    return Case::Runner::runPetscMain(argc, argv);
}
