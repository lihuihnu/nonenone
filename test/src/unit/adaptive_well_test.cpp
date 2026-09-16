/**
 * @file adaptive_well_test.cpp
 * @brief 单元测试：验证 `adaptive_well` 的核心语义、边界条件和回归行为。
 */
#include <adaptive_timestep/well/control_cycle.hpp>

#include <indices/indices.hpp>
#include <well/well.hpp>

#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

using Config = MPMC::CompositionalModelConfig<6, true, true, false, false, false>;
using Indices = MPMC::ScalarIndices<Config>;
using Well = MPMC::WellSpecification<Indices, int>;

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

} // namespace

int main()
{
    Well well{
        1,
        "P1",
        MPMC::WellType::Producer,
        MPMC::WellControl::OilRate,
        100.0,
        0,
        {{0, 1.0}}};

    well.setMinimumBhp(9.0e6);

    std::vector<Well> wells{well};
    std::vector<MPMC::WellState<Indices>> sourceStates(1);
    sourceStates[0].bottomHolePressure = 8.0e6;
    sourceStates[0].surfacePhaseRate[static_cast<std::size_t>(Indices::Phase::liquid)] = -100.0;

    int syncCount = 0;

    auto provider = [&](std::vector<MPMC::WellState<Indices>> &states) {
        states = sourceStates;
    };

    auto sync = [&](const std::vector<Well> &) {
        ++syncCount;
    };

    MPMC::WellControlCycle<Indices, Well, decltype(provider), decltype(sync)> cycle(
        wells,
        provider,
        sync);

    cycle.backup();
    const bool changed = cycle.update();

    require(changed, "Producer minimum-BHP constraint did not switch control.");
    require(wells[0].control == MPMC::WellControl::Bhp, "Control did not switch to BHP.");
    require(wells[0].target == 9.0e6, "BHP target mismatch.");
    require(syncCount == 1, "Control switch did not synchronize wells.");

    cycle.restore();

    require(wells[0].control == MPMC::WellControl::OilRate, "Control snapshot was not restored.");
    require(wells[0].target == 100.0, "Target snapshot was not restored.");
    require(wells[0].primaryControl == MPMC::WellControl::OilRate, "Primary control snapshot changed.");
    require(wells[0].primaryTarget == 100.0, "Primary target snapshot changed.");
    require(wells[0].primaryControlInitialized, "Primary-control initialization flag changed.");
    require(syncCount == 2, "Restore did not synchronize wells.");

    std::cout << "AdaptiveTimeStepper Well control cycle: ALL PASS\n";
    return 0;
}
