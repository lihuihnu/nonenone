/**
 * @file lab_flow_control_contract_test.cpp
 * @brief 回归：实验流量尺度、PVI 与 rate+max-BHP / fixed-BHP 控制契约。
 */
#include "../../../case/scw_kerogen_common/lab_flow_scaling.hpp"

#include <well/control.hpp>
#include <well/specification.hpp>
#include <indices/indices.hpp>
#include <indices/model_config.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
using Config = MPMC::CompositionalModelConfig<
    2, true, true, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Well = MPMC::WellSpecification<Indices, int>;

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
        // Test-only round numbers verify formulas; they are not apparatus values.
        const ScwKerogenLab::FlowScale scale{
            1.0e-4,  // PV [m3]
            1.0e-3,  // area [m2]
            0.20,    // length [m]
            0.25,    // phi
            1.0e-13, // k [m2]
            1.0e-4   // mu [Pa s]
        };
        constexpr double pviRate = 1.0e-4;
        const double q = scale.injectionRateFromPviRate(pviRate);
        require(std::abs(q - 1.0e-8) < 1.0e-20,
                "Q must equal PVI rate times measured PV");
        require(std::abs(scale.nominalResidenceTime(q) - 1.0e4) < 1.0e-8,
                "residence time must equal PV/Q");
        require(std::abs(scale.darcyVelocity(q) - 1.0e-5) < 1.0e-16,
                "Darcy velocity must equal Q/A");
        require(std::abs(scale.poreVelocity(q) - 4.0e-5) < 1.0e-16,
                "pore velocity must equal Q/(phi A)");
        require(std::abs(scale.singlePhaseDarcyPressureDrop(q) - 2000.0) < 1.0e-8,
                "Darcy pressure-drop scale must follow mu L Q/(k A)");
        require(std::abs(ScwKerogenLab::actualPvi(5.0e-5, 1.0e-4) - 0.5) < 1.0e-15,
                "actual PVI must use cumulative actual injected reservoir volume");

        // Injector: primary rate control with maximum-BHP safety switching.
        Well injector(
            0, "INJ", MPMC::WellType::Injector,
            MPMC::WellControl::ReservoirTotalRate,
            q, 0, {});
        injector.setMaximumBhp(29.0e6);
        injector.primaryControl = MPMC::WellControl::ReservoirTotalRate;
        injector.primaryTarget = q;
        injector.primaryControlInitialized = true;

        MPMC::WellState<Indices> injectorState;
        injectorState.bottomHolePressure = 29.2e6;
        injectorState.reservoirPhaseRate[0] = q;
        const auto injUpdate = injector.updateControl(injectorState);
        require(injUpdate.changed,
                "rate injector must switch when maximum BHP is exceeded");
        require(injector.control == MPMC::WellControl::Bhp,
                "injector must switch to BHP control");
        require(std::abs(injector.target - 29.0e6) < 1.0e-6,
                "injector BHP target must equal configured maximum");

        // Producer: fixed BHP is not coupled to injector rate.
        Well producer(
            1, "PROD", MPMC::WellType::Producer,
            MPMC::WellControl::Bhp,
            27.9e6, 1, {});
        producer.primaryControl = MPMC::WellControl::Bhp;
        producer.primaryTarget = 27.9e6;
        producer.primaryControlInitialized = true;

        MPMC::WellState<Indices> producerState;
        producerState.bottomHolePressure = 27.9e6;
        producerState.reservoirPhaseRate[0] = -0.7 * q;
        const auto prodUpdate = producer.updateControl(producerState);
        require(!prodUpdate.changed,
                "fixed-BHP producer must not switch merely to match injector rate");
        require(producer.control == MPMC::WellControl::Bhp,
                "producer must remain fixed BHP");

        std::cout << "Lab flow-control contract: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Lab flow-control contract: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
