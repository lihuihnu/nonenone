/**
 * @file producer_composition_test.cpp
 * @brief Producer composition ledger: mass fractions, cumulative production, RF and E_L/H.
 */
#include <output/metrics/producer_composition.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}
}

int main()
{
    try
    {
        using Ledger = MPMC::ProducerCompositionLedger<5>;
        using Array = Ledger::Array;

        // H2O, Gasoline, Diesel, Middle, Heavy.
        const Array initial{
            10.0, 4.0, 6.0, 10.0, 20.0};

        Ledger ledger;
        ledger.initialize(initial, 0.0);

        // Accepted end-step production rates [kg/s].
        const Array rate1{
            0.10, 0.20, 0.30, 0.20, 0.20};
        ledger.accept(10.0, rate1);

        const std::vector<int> light{1, 2};
        const std::vector<int> heavy{4};
        const auto s1 = ledger.snapshot(rate1, light, heavy);

        require(std::abs(s1.instantaneousMassFraction[0] - 0.10) < 1.0e-14,
                "instantaneous mass fraction mismatch");
        require(std::abs(s1.cumulativeProducedKg[1] - 2.0) < 1.0e-14,
                "cumulative produced mass mismatch");
        require(std::abs(s1.recoveryFraction[1] - 0.5) < 1.0e-14,
                "component recovery fraction mismatch");

        // Produced light/heavy = (0.2+0.3)/0.2 = 2.5.
        // Initial light/heavy = (4+6)/20 = 0.5.
        // Enrichment = 5.
        require(std::abs(s1.instantaneousLightHeavyEnrichment - 5.0) < 1.0e-14,
                "instantaneous light/heavy enrichment mismatch");
        require(std::abs(s1.cumulativeLightHeavyEnrichment - 5.0) < 1.0e-14,
                "cumulative light/heavy enrichment mismatch");

        const Array rate2{
            0.10, 0.10, 0.10, 0.20, 0.50};
        ledger.accept(20.0, rate2);
        const auto s2 = ledger.snapshot(rate2, light, heavy);

        require(std::abs(s2.cumulativeProducedKg[4] - 7.0) < 1.0e-14,
                "second accepted-step cumulative mass mismatch");
        require(std::abs(s2.recoveryFraction[4] - 0.35) < 1.0e-14,
                "heavy recovery fraction mismatch");
        require(s2.instantaneousLightHeavyEnrichment <
                    s1.instantaneousLightHeavyEnrichment,
                "heavier second-step effluent should reduce instantaneous enrichment");

        std::cout << "Producer composition ledger: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Producer composition ledger: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
