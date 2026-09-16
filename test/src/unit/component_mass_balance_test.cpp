/**
 * @file component_mass_balance_test.cpp
 * @brief 单元测试：验证 `component_mass_balance` 的核心语义、边界条件和回归行为。
 */
#include <output/metrics/component_mass_balance.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string &message)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (std::abs(actual - expected) > tolerance * scale)
        throw std::runtime_error(message);
}

void checkBackwardEulerLedger()
{
    using Ledger = MPMC::ComponentMassBalanceLedger<3>;
    using Array = Ledger::Array;

    Ledger ledger;
    require(!ledger.initialized(), "new ledger must be uninitialized");
    ledger.initialize(Array{100.0, 200.0, 300.0}, 0.0);

    // Step 1: dt=10 s, accepted end-state source rates.
    ledger.accept(10.0, Array{2.0, 0.0, 1.0}, Array{0.0, 1.0, 0.5});
    auto snapshot = ledger.snapshot(Array{120.0, 190.0, 305.0});
    for (std::size_t c = 0; c < 3; ++c)
    {
        near(snapshot.error[c], 0.0, 1.0e-14, "step-1 discrete balance");
        near(snapshot.relativeError[c], 0.0, 1.0e-14, "step-1 relative balance");
    }

    // Step 2: dt=5 s, different end-state rates.
    ledger.accept(15.0, Array{0.0, 3.0, 0.0}, Array{4.0, 0.0, 2.0});
    snapshot = ledger.snapshot(Array{100.0, 205.0, 295.0});
    near(snapshot.cumulativeInjected[0], 20.0, 1.0e-14, "CO2 cumulative injection");
    near(snapshot.cumulativeProduced[0], 20.0, 1.0e-14, "CO2 cumulative production");
    near(snapshot.cumulativeInjected[1], 15.0, 1.0e-14, "component-1 cumulative injection");
    near(snapshot.cumulativeProduced[1], 10.0, 1.0e-14, "component-1 cumulative production");
    near(snapshot.cumulativeInjected[2], 10.0, 1.0e-14, "component-2 cumulative injection");
    near(snapshot.cumulativeProduced[2], 15.0, 1.0e-14, "component-2 cumulative production");
    for (std::size_t c = 0; c < 3; ++c)
        near(snapshot.error[c], 0.0, 1.0e-14, "step-2 discrete balance");

    // The reported error must preserve sign and use a finite robust scale.
    snapshot = ledger.snapshot(Array{101.0, 203.0, 295.5});
    near(snapshot.error[0], 1.0, 1.0e-14, "positive balance error sign");
    near(snapshot.error[1], -2.0, 1.0e-14, "negative balance error sign");
    near(snapshot.error[2], 0.5, 1.0e-14, "third balance error");
    for (double value : snapshot.relativeError)
        require(std::isfinite(value), "relative errors must remain finite");
}

void checkGuards()
{
    using Ledger = MPMC::ComponentMassBalanceLedger<2>;
    using Array = Ledger::Array;

    Ledger ledger;
    bool threw = false;
    try
    {
        (void) ledger.snapshot(Array{});
    }
    catch (const std::logic_error &)
    {
        threw = true;
    }
    require(threw, "snapshot before initialize must fail");

    ledger.initialize(Array{1.0, 2.0}, 3.0);
    threw = false;
    try
    {
        ledger.accept(3.0, Array{}, Array{});
    }
    catch (const std::invalid_argument &)
    {
        threw = true;
    }
    require(threw, "accepted times must increase strictly");

    threw = false;
    try
    {
        ledger.accept(4.0, Array{-1.0, 0.0}, Array{});
    }
    catch (const std::invalid_argument &)
    {
        threw = true;
    }
    require(threw, "negative engineering injection magnitude must fail");
}

} // namespace

int main()
{
    try
    {
        checkBackwardEulerLedger();
        checkGuards();
        std::cout << "Component mass-balance ledger: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Component mass-balance ledger: FAIL: " << error.what() << '\n';
        return 1;
    }
}
