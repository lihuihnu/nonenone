/**
 * @file well_control_test.cpp
 * @brief 单元测试：验证 `well_control` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <well/well.hpp>

#include <iostream>
#include <stdexcept>

namespace
{
using Config = MPMC::CompositionalModelConfig<6, true, true>;
using Indices = MPMC::ScalarIndices<Config>;
using Well = MPMC::WellSpecification<Indices, int>;

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void testInjectorSwitching()
{
    Well well(1, "INJ", MPMC::WellType::Injector,
              MPMC::WellControl::TotalRate, 100.0, 0,
              {{0, 1.0e-12}});
    well.setMaximumBhp(20.0e6);

    MPMC::WellState<Indices> state;
    state.bottomHolePressure = 20.2e6;
    state.surfacePhaseRate[Indices::Phase::water] = 100.0;

    const auto toBhp = well.updateControl(state);
    require(toBhp.changed, "injector must switch to maximum BHP");
    require(toBhp.reason == MPMC::WellControlSwitchReason::MaximumBhp,
            "injector switch reason must be MAX_BHP");
    require(toBhp.previousControl == MPMC::WellControl::TotalRate &&
            toBhp.newControl == MPMC::WellControl::Bhp,
            "injector detailed switch controls mismatch");
    require(well.control == MPMC::WellControl::Bhp,
            "injector active BHP control");
    require(well.target == 20.0e6,
            "injector BHP limit target");

    // At maximum BHP, injection capacity above the requested rate means
    // the rate constraint is stricter and must become active again.
    state.bottomHolePressure = 20.0e6;
    state.surfacePhaseRate.fill(0.0);
    state.surfacePhaseRate[Indices::Phase::water] = 120.0;
    require(well.updateControl(state).changed,
            "injector must return to primary rate control");
    require(well.control == MPMC::WellControl::TotalRate,
            "injector primary rate restored");
    require(well.target == 100.0,
            "injector primary target restored");
}

void testProducerSwitching()
{
    Well well(2, "PROD", MPMC::WellType::Producer,
              MPMC::WellControl::OilRate, 50.0, 0,
              {{0, 1.0e-12}});
    well.setMinimumBhp(10.0e6);
    well.setMaximumWaterRate(20.0);

    bool negativeWaterLimitRejected = false;
    try
    {
        well.setMaximumWaterRate(-1.0);
    }
    catch (const std::invalid_argument &)
    {
        negativeWaterLimitRejected = true;
    }
    require(negativeWaterLimitRejected,
            "negative maximum water rate must reject invalid negative limits");

    MPMC::WellState<Indices> state;
    state.bottomHolePressure = 9.8e6;
    state.surfacePhaseRate[Indices::Phase::liquid] = -50.0;
    state.surfacePhaseRate[Indices::Phase::water] = -10.0;

    const auto toBhp = well.updateControl(state);
    require(toBhp.changed, "producer minimum BHP must have priority");
    require(toBhp.reason == MPMC::WellControlSwitchReason::MinimumBhp,
            "producer switch reason must be MIN_BHP");
    require(well.control == MPMC::WellControl::Bhp,
            "producer BHP safety control");

    // At the minimum BHP, excessive water rate is the active rate constraint.
    state.bottomHolePressure = 10.0e6;
    state.surfacePhaseRate[Indices::Phase::liquid] = -45.0;
    state.surfacePhaseRate[Indices::Phase::water] = -25.0;
    require(well.updateControl(state).changed,
            "producer must switch from BHP to water-rate limit");
    require(well.control == MPMC::WellControl::WaterRate,
            "producer water-rate control selected");
    require(well.target == 20.0,
            "producer water-rate target");
}

void testComponentRateReset()
{
    MPMC::WellState<Indices> state;
    state.componentMassRate[0] = -2.5;
    state.componentMassRate[5] = 1.25;
    state.clearRates();
    for (double rate : state.componentMassRate)
        require(rate == 0.0,
                "clearRates must reset per-component well mass rates");
}

} // namespace

int main()
{
    testInjectorSwitching();
    testProducerSwitching();
    testComponentRateReset();
    std::cout << "Well control validation: ALL PASS\n";
    return 0;
}
