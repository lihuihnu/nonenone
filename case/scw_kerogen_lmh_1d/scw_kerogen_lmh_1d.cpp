/**
 * @file scw_kerogen_lmh_1d.cpp
 * @brief 第二阶段超临界水驱轻中重裂解产物拟组分算例入口。
 */
#include <case/petsc_custom_hooks.hpp>
#include <case/structured_single_eos_reservoir_runner.hpp>

#include "case_config.hpp"
#include "well_config.hpp"

namespace Case
{
struct Definition final
{
    using Config = CaseConfig::Config;
    [[nodiscard]] static constexpr const auto &wells() noexcept
    {
        return WellConfig::wells;
    }
};
using Runner = MPMC::cases::StructuredSingleEosReservoirRunner<Definition>;
using Runtime = typename Runner::Runtime;
} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    return Case::Runner::runPetscMain(argc, argv);
}
