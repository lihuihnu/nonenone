/**
 * @file sw_flash_recovery_test.cpp
 * @brief 单元测试：验证 `sw_flash_recovery` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <indices/model_config.hpp>
#include <natural/compositional_mixture.hpp>
#include <natural/thermo/cubic_eos.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
using Config = MPMC::CompositionalModelConfig<
    4, true, false, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using Indices = MPMC::ScalarIndices<Config>;
using Eos = MPMC::CubicEquationOfState<Indices>;
using Flash = MPMC::CubicThreePhaseFlash<Indices>;
using Composition = std::array<double, 4>;

constexpr std::array<double, 4> tc{647.096, 304.1282, 190.564, 723.0};
constexpr std::array<double, 4> pc{22.064e6, 7.3773e6, 4.5992e6, 1.410e6};
constexpr std::array<double, 4> vc{5.6e-5, 9.4e-5, 9.9e-5, 9.0e-4};
constexpr std::array<double, 4> om{0.3443, 0.22394, 0.01142, 0.742};
constexpr std::array<double, 4> mw{0.01801528, 0.0440095, 0.016043, 0.226441};
constexpr std::array<std::array<double, 4>, 4> kij{{
    {{0.0, 0.1896, 0.485, 0.5}},
    {{0.1896, 0.0, 0.12, 0.09}},
    {{0.485, 0.12, 0.0, 0.0}},
    {{0.5, 0.09, 0.0, 0.0}}
}};

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

Eos makeSw()
{
    MPMC::CompositionalMixture<Indices> mixture(tc, pc, vc, om, mw, kij);
    Eos eos(0.45724,
            0.07780,
            std::move(mixture),
            5.0,
            2.414213562373095,
            -0.414213562373095,
            1e-30);

    Eos::SoreideWhitsonOptions sw;
    sw.waterComponent = 0;
    sw.salinityMolality = 0.0;
    sw.aqueousWaterBip[0] = [](double, double) { return 0.0; };
    sw.aqueousWaterBip[1] = [](double temperature, double salinity) {
        return MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
            temperature, tc[1], salinity);
    };
    sw.aqueousWaterBip[2] = [](double temperature, double salinity) {
        return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
            temperature, tc[2], om[2], salinity);
    };
    sw.aqueousWaterBip[3] = [](double, double) { return 0.5; };
    eos.configureSoreideWhitson(std::move(sw));
    return eos;
}

void checkMassBalance(const MPMC::ThreePhaseFlashResult<Indices>& result,
                      const Composition& overall,
                      const std::string& name)
{
    for (std::size_t i = 0; i < overall.size(); ++i) {
        double reconstructed = 0.0;
        for (int phase = 0; phase < 3; ++phase) {
            if (!result.presence.contains(static_cast<MPMC::CompositionalPhase>(phase))) {
                continue;
            }
            reconstructed += result.phaseMoleFraction[phase] *
                             result.composition[phase][i];
        }
        require(std::abs(reconstructed - overall[i]) < 2e-12,
                name + ": component mass balance");
    }
}

void checkFugacityClosure(const Eos& eos,
                          double pressure,
                          double temperature,
                          const MPMC::ThreePhaseFlashResult<Indices>& result,
                          const std::string& name)
{
    std::array<Eos::PhaseResult<double>, 3> phaseResult{};
    std::array<bool, 3> active{false, false, false};
    const std::array<MPMC::CompositionalPhase, 3> roles{
        MPMC::CompositionalPhase::Oil,
        MPMC::CompositionalPhase::Gas,
        MPMC::CompositionalPhase::Water};

    for (int phase = 0; phase < 3; ++phase) {
        if (!result.presence.contains(static_cast<MPMC::CompositionalPhase>(phase))) {
            continue;
        }
        phaseResult[phase] = eos.phaseResult(
            pressure, temperature, result.composition[phase], roles[phase]);
        active[phase] = true;
    }

    for (std::size_t i = 0; i < 4; ++i) {
        double lo = 1e300;
        double hi = -1e300;
        int count = 0;
        for (int phase = 0; phase < 3; ++phase) {
            if (!active[phase] || result.composition[phase][i] <= 1e-20) {
                continue;
            }
            const double logF = std::log(std::max(
                phaseResult[phase].fugacity[i], 1e-300));
            lo = std::min(lo, logF);
            hi = std::max(hi, logF);
            ++count;
        }
        if (count > 1) {
            require(hi - lo < 2e-9, name + ": fugacity closure");
        }
    }
}

void checkStableReducedSet(const Eos& eos)
{
    constexpr Composition z{0.10, 0.25, 0.50, 0.15};
    constexpr double temperature = 280.0;
    constexpr double pressureBar = 262.567915178188;
    const double pressure = pressureBar * 1e5;

    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    Flash flash(eos, options);
    const auto result = flash.flash(pressure, temperature, z);

    require(result.converged, "SW stable-reduced fallback: convergence");
    require(result.presence.bits() ==
                (MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit),
            "SW stable-reduced fallback: expected O+W");
    require(result.phaseMoleFraction[1] == 0.0,
            "SW stable-reduced fallback: gas must be absent");
    require(result.composition[2][0] > result.composition[0][0],
            "SW stable-reduced fallback: water identity");

    checkMassBalance(result, z, "SW stable-reduced fallback");
    checkFugacityClosure(eos,
                         pressure,
                         temperature,
                         result,
                         "SW stable-reduced fallback");
}

void checkAllocationThreePhase(const Eos& eos)
{
    constexpr Composition z{0.10, 0.25, 0.50, 0.15};
    constexpr double temperature = 404.0;
    constexpr double pressureBar = 16.203947518373063;
    const double pressure = pressureBar * 1e5;

    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    Flash flash(eos, options);
    const auto result = flash.flash(pressure, temperature, z);

    require(result.converged, "SW allocation fallback: convergence");
    // The raw allocation basin has coincident O/G compositions (maximum
    // component difference below 2e-11).  They are one thermodynamic phase,
    // so the public result must combine their mole fractions exactly instead
    // of exposing a numerical three-phase duplicate.
    require(result.presence.bits() ==
                (MPMC::PhasePresence::oilBit | MPMC::PhasePresence::waterBit),
            "SW allocation fallback: coincident O/G must collapse to O+W");
    require(result.phaseMoleFraction[0] > 0.15 &&
                result.phaseMoleFraction[1] == 0.0 &&
                result.phaseMoleFraction[2] > 0.80,
            "SW allocation fallback: collapsed physical phase fractions");
    require(result.composition[2][0] > result.composition[0][0] &&
                result.composition[2][0] > result.composition[1][0],
            "SW allocation fallback: water identity");

    checkMassBalance(result, z, "SW allocation fallback");
    checkFugacityClosure(eos,
                         pressure,
                         temperature,
                         result,
                         "SW allocation fallback");
}

void checkPressureContinuationThreePhase(const Eos& eos)
{
    // v41's final P-T-map failure: the physical O+G+W solution exists and is
    // smooth from the lower-pressure branch, but default Wilson/water K seeds
    // miss its attraction basin at this pressure.  The fail-only pressure
    // homotopy added in v42 must recover the equilibrium without changing EOS
    // parameters, tolerances, or public phase labels.
    constexpr Composition z{0.10, 0.25, 0.50, 0.15};
    constexpr double temperature = 300.0;
    constexpr double pressureBar = 206.59796623849141;
    const double pressure = pressureBar * 1e5;

    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = 0;
    Flash flash(eos, options);
    const auto result = flash.flash(pressure, temperature, z);

    require(result.converged, "SW pressure-continuation fallback: convergence");
    require(result.presence.bits() == MPMC::PhasePresence::allBits,
            "SW pressure-continuation fallback: expected O+G+W");
    require(result.phaseMoleFraction[0] > 0.70 &&
                result.phaseMoleFraction[1] > 0.15 &&
                result.phaseMoleFraction[2] > 0.09,
            "SW pressure-continuation fallback: physical phase fractions");
    require(result.composition[2][0] > result.composition[0][0] &&
                result.composition[2][0] > result.composition[1][0],
            "SW pressure-continuation fallback: water identity");

    checkMassBalance(result, z, "SW pressure-continuation fallback");
    checkFugacityClosure(eos,
                         pressure,
                         temperature,
                         result,
                         "SW pressure-continuation fallback");
}
} // namespace

int main()
{
    try {
        const auto eos = makeSw();
        checkStableReducedSet(eos);
        checkAllocationThreePhase(eos);
        checkPressureContinuationThreePhase(eos);
        std::cout << "SW fail-only flash recovery: ALL PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SW fail-only flash recovery: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
