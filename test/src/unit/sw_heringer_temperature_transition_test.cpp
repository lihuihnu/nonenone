/**
 * @file sw_heringer_temperature_transition_test.cpp
 * @brief 回归验证 Heringer BSB 算例中 SW 升温扫描的相态拓扑。
 */
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
using Config = MPMC::CompositionalModelConfig<
    8, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Composition = std::array<double, 8>;

constexpr Composition tc{
    647.30, 304.20, 160.00, 344.22, 463.22, 605.78, 751.00, 942.50};
constexpr Composition pc{
    220.48e5, 73.76e5, 46.00e5, 45.00e5,
    34.00e5, 21.75e5, 16.54e5, 16.42e5};
constexpr Composition vc{
    5.60e-5, 9.43e-5, 9.93e-5, 1.81e-4,
    3.06e-4, 5.99e-4, 1.13e-3, 2.09e-3};
constexpr Composition omega{
    0.344, 0.225, 0.008, 0.131, 0.240, 0.618, 0.957, 1.268};
constexpr Composition mw{
    0.01801528, 0.044010, 0.016040, 0.037200,
    0.069500, 0.140960, 0.280990, 0.519620};
constexpr Composition feed{
    0.750000, 0.105055, 0.012915, 0.022545,
    0.025065, 0.049560, 0.024165, 0.010695};
constexpr std::array<std::array<double, 8>, 8> kij{{
    {{0.0000, 0.1896, 0.4850, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000}},
    {{0.1896, 0.0000, 0.0550, 0.0550, 0.0550, 0.1050, 0.1050, 0.1050}},
    {{0.4850, 0.0550, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.0550, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.0550, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.1050, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.1050, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}},
    {{0.5000, 0.1050, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}}
}};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

Eos makeStandardSw()
{
    Eos eos(
        0.45724,
        0.07780,
        MPMC::CompositionalMixture<Indices>(tc, pc, vc, omega, mw, kij),
        5.0,
        2.414213562373095,
        -0.414213562373095,
        1.0e-30);
    Eos::SoreideWhitsonOptions sw;
    sw.waterComponent = 0;
    sw.salinityMolality = 0.0;
    sw.aqueousWaterBip[0] = [](double, double) { return 0.0; };
    // Empty non-water callbacks retain Heringer Table B4 BIPs.  This is the
    // standard, non-fitted SW case used by the publication scan.
    eos.configureSoreideWhitson(std::move(sw));
    return eos;
}

void checkMassBalance(const Flash::Result& result, const std::string& name)
{
    for (std::size_t component = 0; component < feed.size(); ++component) {
        double reconstructed = 0.0;
        for (std::size_t phase = 0; phase < 3; ++phase) {
            reconstructed += result.phaseMoleFraction[phase] *
                             result.composition[phase][component];
        }
        require(std::abs(reconstructed - feed[component]) < 3.0e-12,
                name + ": component mass balance");
    }
}

void checkState(const Flash& flash,
                double pressureMpa,
                double temperature,
                unsigned expectedPresence,
                double minimumOilFraction,
                const std::string& name)
{
    const double pressure = pressureMpa * 1.0e6;
    const auto result = flash.flash(pressure, temperature, feed);
    require(result.converged, name + ": convergence");
    require(result.presence.bits() == expectedPresence,
            name + ": phase topology");
    require(result.phaseMoleFraction[0] >= minimumOilFraction,
            name + ": oil fraction must not collapse to metastable water");
    checkMassBalance(result, name);
}
} // namespace

int main()
{
    try {
        const Eos eos = makeStandardSw();
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = 0;
        options.maximumIterations = 240;
        options.maximumStabilityIterations = 160;
        options.fugacityTolerance = 1.0e-10;
        const Flash flash(eos, options);

        constexpr unsigned ow =
            MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit;
        checkState(flash, 28.0, 460.0, ow, 0.25, "SW at 460 K");
        checkState(flash, 28.0, 500.0, ow, 0.26, "SW at 500 K");
        checkState(flash, 28.0, 530.0, ow, 0.28, "SW at 530 K");
        checkState(flash, 28.0, 550.0, MPMC::PhasePresence::allBits, 0.25,
                   "SW at 550 K");
        checkState(flash, 28.0, 630.0, ow, 0.04, "SW at 630 K");
        checkState(flash, 28.0, 690.0, MPMC::PhasePresence::waterBit, 0.0,
                   "SW at 690 K");
        checkState(flash, 6.0, 300.0, MPMC::PhasePresence::allBits, 0.20,
                   "SW P-T topology at 300 K and 6 MPa");
        checkState(flash, 33.5, 600.0, MPMC::PhasePresence::allBits, 0.10,
                   "SW P-T topology at 600 K and 33.5 MPa");

        std::cout << "SW Heringer temperature transition: ALL PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SW Heringer temperature transition: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
