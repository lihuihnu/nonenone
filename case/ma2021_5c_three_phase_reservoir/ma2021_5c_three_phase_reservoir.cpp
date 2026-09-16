/**
 * @file ma2021_5c_three_phase_reservoir.cpp
 * @brief Ma 2021 五组分油气水三相文献算例入口。
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
