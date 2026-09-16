/**
 * @file three_eos_3d_case_test.cpp
 * @brief 单元测试：验证三维 PR/SW/CPA 对比算例的共同物理配置和三套初始 flash。
 */
#include "../../../case/three_eos_3d_compare/case_config.hpp"
#include "../../../case/three_eos_3d_compare/well_config.hpp"

#include <case/well_factory.hpp>
#include <case/case_validation.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <natural/state/three_phase_equilibrium.hpp>
#include <natural/numerics.hpp>
#include <natural/petsc/phase_state_codec.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells,
    CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption,
    CaseConfig::Model::enableLandTrapping,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;

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

template <class FactoryConfig>
auto flashInitial()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, FactoryConfig>();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = CaseConfig::CommonFluid::waterComponent;
    const MPMC::CubicThreePhaseFlash<Indices> flash(fluid.eos, options);
    return flash.flash(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);
}

void checkPhaseSuppressionCodec()
{
    MPMC::PhaseStateData<Indices> state;
    state.phasePresence = MPMC::PhasePresence::gasOnly();
    state.phasePresence.add(MPMC::CompositionalPhase::Water);
    state.phaseSuppression = MPMC::PhasePresence(std::uint8_t{0});
    state.phaseSuppression.add(MPMC::CompositionalPhase::Oil);

    const auto encoded = MPMC::PetscPhaseStateCodec<Indices>::encode(state);
    const auto decoded = MPMC::PetscPhaseStateCodec<Indices>::decode(encoded);
    require(decoded.phasePresence.bits() == state.phasePresence.bits(),
            "phase-state codec must preserve public phase presence");
    require(decoded.phaseSuppression.bits() == state.phaseSuppression.bits(),
            "phase-state codec must preserve hysteresis suppression memory");

    require(encoded[static_cast<std::size_t>(Indices::PhaseState::flag)] ==
                MPMC::encodePhasePresence(state.phasePresence),
            "public phase_presence_mask must remain the historical 1..7 encoding");
    require(encoded[static_cast<std::size_t>(Indices::PhaseState::phaseSuppressionFlag)] ==
                static_cast<double>(state.phaseSuppression.bits()),
            "hysteresis suppression must use its own phase-state column");
}

void checkConfiguration()
{
    static_assert(Indices::fullyCompositionalThreePhase);
    static_assert(!Indices::hasAqueousCO2Dissolution);
    static_assert(!Indices::hasAdsorption);
    static_assert(!Indices::hasLandTrapping);
    static_assert(CaseConfig::Grid::nz > 1);
    static_assert(CaseConfig::Model::numberOfComponents == 5);
    static_assert(CaseConfig::Time::adaptive,
                  "three-EOS reservoir comparison requires transactional adaptive substeps");
    static_assert(CaseConfig::Time::printAdaptiveSteps,
                  "three-EOS diagnostic benchmark must print ACCEPT/REJECT/RETRY events");
    static_assert(CaseConfig::Output::every == 1,
                  "diagnostic benchmark must print every common output target");
    static_assert(CaseConfig::Output::printWells);
    static_assert(CaseConfig::Output::printWellPhaseDetails);
    static_assert(CaseConfig::Output::printInventory);
    static_assert(CaseConfig::Output::printComponentMassBalance);
    static_assert(CaseConfig::Output::printNewtonIterations);
    static_assert(CaseConfig::Output::printWellControlSwitches);
    static_assert(CaseConfig::Output::printReservoirDiagnostics);
    static_assert(CaseConfig::Output::printFinalSummary);
    static_assert(CaseConfig::PrFluid::thermodynamicModel ==
                  MPMC::CubicThermodynamicModel::PengRobinson);
    static_assert(CaseConfig::SwFluid::thermodynamicModel ==
                  MPMC::CubicThermodynamicModel::SoreideWhitson);
    static_assert(CaseConfig::CpaFluid::thermodynamicModel ==
                  MPMC::CubicThermodynamicModel::CubicPlusAssociation);

    MPMC::cases::validateCaseConfig<Indices, CaseConfig::Config>();
    near(CaseConfig::InitialState::pressure / CaseConfig::bar, 60.0, 1.0e-13,
         "common initial pressure");
    near(CaseConfig::InitialState::temperature, 305.0, 1.0e-13,
         "common initial temperature");

    double sum = 0.0;
    for (double z : CaseConfig::InitialState::overallComposition)
        sum += z;
    near(sum, 1.0, 1.0e-13, "overall composition normalization");

    // 公开数据锚点用于防止算例参数在后续维护中偏离已记录的 CoolProp/SW/CPA 来源。
    near(CaseConfig::CommonFluid::criticalTemperature[0], 647.096, 1.0e-12,
         "CoolProp water critical temperature");
    near(CaseConfig::CommonFluid::criticalTemperature[4], 425.125, 1.0e-12,
         "CoolProp n-butane critical temperature");
    near(CaseConfig::SwFluid::soreideWhitsonSalinityMolality, 0.0, 1.0e-13,
         "SW fresh-water salinity");
    near(CaseConfig::CpaFluid::cpaB[0], 14.52e-6, 1.0e-13,
         "CPA water co-volume");
    near(CaseConfig::CpaFluid::cpaGamma[4], 2193.08, 1.0e-12,
         "CPA n-butane Gamma");
    near(CaseConfig::CpaFluid::cpaAssociationVolume[0], 0.0692, 1.0e-13,
         "CPA 4C-water association volume");

    require(WellConfig::wells.size() == 2, "comparison case requires one injector and one producer");
    near(WellConfig::wells[0].target * CaseConfig::secondsPerDay, 300000.0, 1.0e-13,
         "CO2 injector surface/reference rate");
    require(CaseConfig::Rock::kx(3, 4, 0) > CaseConfig::Rock::kx(3, 20, 0),
            "injector should lie inside the high-permeability channel");
    require(CaseConfig::Rock::kz(18, CaseConfig::Rock::channelCenter(18), 3) >
                CaseConfig::Rock::kz(5, CaseConfig::Rock::channelCenter(5), 3),
            "baffle must contain a localized vertical-flow window");
}


void checkPrNearPhaseBoundaryHysteresis()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::PrFactoryConfig>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    // Regression state reconstructed from the v54 8-rank PR near-well failure:
    // the unrestricted P-T-z flash contains only an O(1e-6) oil-rich phase.
    // A GW restricted flash is only marginally unstable, so the simulator
    // active set must cross the phase boundary instead of asymptotically
    // shrinking dt while retaining an almost-zero oil phase.
    const std::array<double, Indices::numComponents> z{
        0.42492349, 0.52532524, 0.00465313, 0.00887512, 0.03622302};
    const double pressure = 61.0917 * CaseConfig::bar;
    const auto flash = equilibrium.flashPTZ(
        pressure, CaseConfig::InitialState::temperature, z);

    require(flash.converged, "near-boundary PR flash must converge");
    require(flash.presence.count() == 3,
            "near-boundary unrestricted PR flash should still report O/G/W");
    require(flash.saturation[0] > 0.0 &&
                flash.saturation[0] < MPMC::NaturalNumerics::phaseDisappearanceSaturation,
            "near-boundary oil saturation must lie inside disappearance hysteresis");

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = pressure;
    primary[Indices::Primary::wellPressure] = pressure;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, flash);
    equilibrium.updatePhaseState(primary, phaseState);

    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "marginal O(1e-6) oil phase must disappear from the simulator active set");
    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Gas) &&
                phaseState.phasePresence.contains(MPMC::CompositionalPhase::Water),
            "near-boundary state must reduce to G+W");
    near(primary[Indices::Primary::liquidSaturation], 0.0, 0.0,
         "inactive oil saturation must be exactly zero");
}

void checkSwNearPhaseBoundaryHysteresis()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::SwFactoryConfig>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    // Regression state reconstructed from the v55 target-HPC SW failure at
    // t=8.65179857 day.  The oil-rich phase carries only O(1e-6) saturation
    // and O(1e-8) mobility, while the time-step controller reaches its minimum
    // dt before the old 1e-6 disappearance trigger can be crossed.
    constexpr double pressure = 61.0907274 * CaseConfig::bar;
    constexpr std::array<double, 3> saturation{
        1.63784632e-6,
        0.929539171,
        0.0704591913};
    constexpr std::array<std::array<double, Indices::numComponents>, 3> composition{{
        {{0.00215254, 0.82912252, 0.00394291, 0.01458792, 0.15019411}},
        {{0.00151687, 0.91213046, 0.00807785, 0.01539587, 0.06287895}},
        {{0.982817561, 0.0171256784, 1.24791813e-5, 1.96386947e-5, 2.46430195e-5}}
    }};

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = pressure;
    primary[Indices::Primary::wellPressure] = pressure;
    primary[Indices::Primary::liquidSaturation] = saturation[0];
    primary[Indices::Primary::vaporSaturation] = saturation[1];
    primary[Indices::Primary::waterSaturation] = saturation[2];

    const auto writeComposition = [&](const auto &indices, const auto &values) {
        for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
            primary[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])] =
                values[static_cast<std::size_t>(i)];
    };
    writeComposition(Indices::Primary::liquidComposition, composition[0]);
    writeComposition(Indices::Primary::vaporComposition, composition[1]);
    writeComposition(Indices::Primary::waterComposition, composition[2]);

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence::all();
    equilibrium.updateSecondary(primary, phaseState);

    require(primary[Indices::Primary::liquidSaturation] > 0.0 &&
                primary[Indices::Primary::liquidSaturation] <
                    MPMC::NaturalNumerics::phaseDisappearanceSaturation,
            "v56 SW trace oil must enter the reduced-phase probe band");

    equilibrium.updatePhaseState(primary, phaseState);

    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "v56 SW trace oil must disappear after guarded reduced-phase test");
    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Gas) &&
                phaseState.phasePresence.contains(MPMC::CompositionalPhase::Water),
            "v56 SW boundary state must reduce to G+W");
    near(primary[Indices::Primary::liquidSaturation], 0.0, 0.0,
         "inactive SW oil saturation must be exactly zero");
}


void checkConsistentPhaseBoundaryHysteresis()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::SwFactoryConfig>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    // v57 SW 目标超算在第二个相边界反复出现 So≈1.1e-5：旧实现先用
    // fail-only 1e-5 stability deadband 删除油相，随后普通 updateState 又用
    // 1e-6 deadband 立即把油相重新生成，导致同一 active set 来回切换并把 dt
    // 压到 minimumDt。现在所有 Newton/update/retry 共用唯一一套滞回规则。
    constexpr double pressure = 61.0907255 * CaseConfig::bar;
    constexpr std::array<double, 3> saturation{
        1.10e-5,
        0.9295298,
        1.0 - 0.9295298 - 1.10e-5};
    constexpr std::array<std::array<double, Indices::numComponents>, 3> composition{{
        {{0.00215254, 0.82912252, 0.00394291, 0.01458792, 0.15019411}},
        {{0.00151687, 0.91213046, 0.00807785, 0.01539587, 0.06287895}},
        {{0.982817561, 0.0171256784, 1.24791813e-5, 1.96386947e-5, 2.46430195e-5}}
    }};

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = pressure;
    primary[Indices::Primary::wellPressure] = pressure;
    primary[Indices::Primary::liquidSaturation] = saturation[0];
    primary[Indices::Primary::vaporSaturation] = saturation[1];
    primary[Indices::Primary::waterSaturation] = saturation[2];
    const auto writeComposition = [&](const auto &indices, const auto &values) {
        for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
            primary[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])] =
                values[static_cast<std::size_t>(i)];
    };
    writeComposition(Indices::Primary::liquidComposition, composition[0]);
    writeComposition(Indices::Primary::vaporComposition, composition[1]);
    writeComposition(Indices::Primary::waterComposition, composition[2]);

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence::all();
    equilibrium.updateSecondary(primary, phaseState);

    require(primary[Indices::Primary::liquidSaturation] <
                MPMC::NaturalNumerics::phaseDisappearanceSaturation,
            "v57 SW plateau oil must enter the unified disappearance probe band");

    equilibrium.updatePhaseState(primary, phaseState);
    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "unified active set must reduce the v57 SW plateau state to G+W");
    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Gas) &&
                phaseState.phasePresence.contains(MPMC::CompositionalPhase::Water),
            "unified active set must retain G+W");
    near(primary[Indices::Primary::liquidSaturation], 0.0, 0.0,
         "inactive oil saturation must be exactly zero");

    // 最关键回归：相已经消失后再次执行普通 updateState，不允许因为换回另一套
    // appearance margin 而立刻重新出现。这里直接重复调用来模拟下一次 Newton。
    equilibrium.updatePhaseState(primary, phaseState);
    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "phase-boundary hysteresis must be idempotent across consecutive Newton updates");
    near(primary[Indices::Primary::liquidSaturation], 0.0, 0.0,
         "repeated active-set update must keep inactive oil saturation at zero");
}

template <class FactoryConfig>
void checkLateBoundaryState(
    const char *label,
    double pressureBar,
    const std::array<double, 3> &saturation,
    const std::array<std::array<double, Indices::numComponents>, 3> &composition)
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, FactoryConfig>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = pressureBar * CaseConfig::bar;
    primary[Indices::Primary::wellPressure] = pressureBar * CaseConfig::bar;
    primary[Indices::Primary::liquidSaturation] = saturation[0];
    primary[Indices::Primary::vaporSaturation] = saturation[1];
    primary[Indices::Primary::waterSaturation] = saturation[2];

    const auto writeComposition = [&](const auto &indices, const auto &values) {
        for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
            primary[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])] =
                values[static_cast<std::size_t>(i)];
    };
    writeComposition(Indices::Primary::liquidComposition, composition[0]);
    writeComposition(Indices::Primary::vaporComposition, composition[1]);
    writeComposition(Indices::Primary::waterComposition, composition[2]);

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence::all();
    equilibrium.updateSecondary(primary, phaseState);

    require(primary[Indices::Primary::liquidSaturation] <
                MPMC::NaturalNumerics::phaseDisappearanceSaturation,
            std::string(label) + " oil saturation must enter unified probe band");

    equilibrium.updatePhaseState(primary, phaseState);
    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            std::string(label) + " must reduce from OGW to GW");
    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Gas) &&
                phaseState.phasePresence.contains(MPMC::CompositionalPhase::Water),
            std::string(label) + " reduced state must be GW");

    // 再做一次普通 Newton state update，专门防止 v57 的“删相后立即再生”回归。
    equilibrium.updatePhaseState(primary, phaseState);
    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            std::string(label) + " GW state must remain stable on repeated update");
}

void checkLatePhaseBoundaryPlateausAcrossEos()
{
    // v57 SW: t≈29.40495 day, cell 183。日志给的是质量分数；这里已按算例
    // molar mass 转回摩尔分数并归一化。
    checkLateBoundaryState<CaseConfig::SwFactoryConfig>(
        "late SW plateau",
        60.8573193,
        {{1.12379813e-5, 0.928707307, 0.0712814551}},
        {{{{0.00214849135469, 0.827066798686, 0.00369369301805, 0.0142636880279, 0.152827328914}},
          {{0.00151582689914, 0.912535784248, 0.00763269414199, 0.0150526101581, 0.0632630845527}},
          {{0.982850464714, 0.0170937878335, 1.17247501635e-5, 1.91554930227e-5, 2.48672090423e-5}}}});

    // v57 CPA: t≈15.64603 day, cell 1011。
    checkLateBoundaryState<CaseConfig::CpaFactoryConfig>(
        "late CPA plateau",
        61.2119752,
        {{1.24332781e-5, 0.941118346, 0.0588692204}},
        {{{{0.000314811660067, 0.845842765216, 0.000998391655232, 0.0067975258748, 0.146046505594}},
          {{0.00072272582824, 0.940153749256, 0.0024644587276, 0.00735180698013, 0.0493072592085}},
          {{0.996355796513, 0.00364383660695, 2.50253195586e-7, 1.09582065508e-7, 7.04483715996e-9}}}});
}


void checkV58CpaMinimumDtPlateauRegression()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::CpaFactoryConfig>();
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    // Reconstructed from job 5945345 at t=10.3435286367 day, cell 147.
    // v58 stopped at minimumDt with So=2.11509217e-5: this is just outside the
    // old 2e-5 probe, while the reduced G+W stability trial sum is only weakly
    // unstable (O(1e-5)).  The v59 numerical deadband must treat this as a
    // trace-phase active-set boundary instead of forcing another dt cut.
    constexpr double pressure = 61.0937428 * CaseConfig::bar;
    constexpr std::array<double, 3> saturation{
        2.11509217e-5,
        0.941339392,
        0.0586394566};
    constexpr std::array<std::array<double, Indices::numComponents>, 3> composition{{
        {{0.000313478202066, 0.842418526326, 0.00122972407829,
          0.00750501899777, 0.148533252396}},
        {{0.000724310348298, 0.938385239068, 0.00304561259234,
          0.00810146029561, 0.0497433776955}},
        {{0.996366841117, 0.00363272290524, 3.08336879942e-7,
          1.20532127501e-7, 7.10906910153e-9}}
    }};

    require(saturation[0] > 2.0e-5,
            "v58 CPA regression state must stay outside the historical 2e-5 probe");
    require(saturation[0] < MPMC::NaturalNumerics::phaseBoundaryProbeSaturation,
            "v59 CPA regression state must enter the thermodynamic probe band");

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = pressure;
    primary[Indices::Primary::wellPressure] = pressure;
    primary[Indices::Primary::liquidSaturation] = saturation[0];
    primary[Indices::Primary::vaporSaturation] = saturation[1];
    primary[Indices::Primary::waterSaturation] = saturation[2];

    const auto writeComposition = [&](const auto &indices, const auto &values) {
        for (int i = 0; i < Indices::numIndependentCompositionsPerPhase; ++i)
            primary[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])] =
                values[static_cast<std::size_t>(i)];
    };
    writeComposition(Indices::Primary::liquidComposition, composition[0]);
    writeComposition(Indices::Primary::vaporComposition, composition[1]);
    writeComposition(Indices::Primary::waterComposition, composition[2]);

    MPMC::PhaseStateData<Indices> phaseState;
    phaseState.phasePresence = MPMC::PhasePresence::all();
    equilibrium.updateSecondary(primary, phaseState);
    equilibrium.updatePhaseState(primary, phaseState);

    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "v58 minimum-dt CPA plateau must reduce to G+W in the v59 active set");
    require(phaseState.phaseSuppression.contains(MPMC::CompositionalPhase::Oil),
            "CPA boundary reduction must retain hysteresis suppression history");
    near(primary[Indices::Primary::liquidSaturation], 0.0, 0.0,
         "v59 CPA boundary reduction must set inactive oil saturation to zero");

    equilibrium.updatePhaseState(primary, phaseState);
    require(!phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "v59 CPA reduced state must remain idempotent on the next Newton update");
}

void checkStrongOilInstabilityReappears()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::PrFactoryConfig>();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = CaseConfig::CommonFluid::waterComponent;
    MPMC::CubicThreePhaseFlash<Indices> flash(fluid.eos, options);
    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);

    auto gw = MPMC::PhasePresence::gasOnly();
    gw.add(MPMC::CompositionalPhase::Water);
    const auto reduced = flash.flashRestricted(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition,
        gw);
    require(reduced.converged && reduced.presence.bits() == gw.bits(),
            "initial P-T-z must admit a restricted G+W equilibrium seed");

    const auto stability = flash.stabilityTest(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition,
        reduced.presence,
        reduced.composition);
    require(stability.missingPhaseUnstable[0] &&
                stability.trialSum[0] >
                    1.0 + MPMC::NaturalNumerics::phaseAppearanceStabilityMargin,
            "strong missing-oil instability must exceed active-set hysteresis");

    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = CaseConfig::InitialState::pressure;
    primary[Indices::Primary::wellPressure] = CaseConfig::InitialState::pressure;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, reduced);
    equilibrium.updatePhaseState(primary, phaseState);

    require(phaseState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "clearly unstable missing oil phase must reappear");
    require(primary[Indices::Primary::liquidSaturation] >
                MPMC::NaturalNumerics::phaseDisappearanceSaturation,
            "reappeared oil phase must receive a physical non-trace saturation");

    // 即使这个缺失油相带有 hysteresis suppression 历史，只要稳定性明确失稳，
    // 仍必须恢复三相并清除 suppression；滞回绝不能变成永久锁相。
    auto suppressedPrimary = primary;
    auto suppressedState = phaseState;
    equilibrium.assignFlashResult(suppressedPrimary, suppressedState, reduced);
    suppressedState.phaseSuppression.add(MPMC::CompositionalPhase::Oil);
    equilibrium.updatePhaseState(suppressedPrimary, suppressedState);
    require(suppressedState.phasePresence.contains(MPMC::CompositionalPhase::Oil),
            "strongly unstable suppressed oil phase must still reappear");
    require(!suppressedState.phaseSuppression.contains(MPMC::CompositionalPhase::Oil),
            "reappeared oil phase must clear hysteresis suppression memory");
}

void checkInitialFlash()
{
    const auto pr = flashInitial<CaseConfig::PrFactoryConfig>();
    const auto sw = flashInitial<CaseConfig::SwFactoryConfig>();
    const auto cpa = flashInitial<CaseConfig::CpaFactoryConfig>();

    require(pr.converged && sw.converged && cpa.converged,
            "all three EOS must converge at the common initial P-T-z state");
    for (const auto *result : {&pr, &sw, &cpa})
    {
        require(result->phaseMoleFraction[0] > 1.0e-4,
                "initial state must contain an oil-rich phase");
        require(result->phaseMoleFraction[1] > 1.0e-4,
                "initial state must contain a gas-rich phase");
        require(result->phaseMoleFraction[2] > 1.0e-4,
                "initial state must contain a water-rich phase");
        near(result->saturation[0] + result->saturation[1] + result->saturation[2],
             1.0, 1.0e-11, "saturation closure");
    }

    // Regression anchors generated by the same production flash/factory path.
    near(pr.saturation[0], 0.5937905132, 2.0e-8, "PR initial So");
    near(pr.saturation[1], 0.3329129536, 2.0e-8, "PR initial Sg");
    near(pr.saturation[2], 0.07329653322, 2.0e-8, "PR initial Sw");

    near(sw.saturation[0], 0.5949988308, 2.0e-8, "SW initial So");
    near(sw.saturation[1], 0.3309409565, 2.0e-8, "SW initial Sg");
    near(sw.saturation[2], 0.07406021268, 2.0e-8, "SW initial Sw");

    near(cpa.saturation[0], 0.5626431812, 2.0e-8, "CPA initial So");
    near(cpa.saturation[1], 0.3778841217, 2.0e-8, "CPA initial Sg");
    near(cpa.saturation[2], 0.05947269707, 2.0e-8, "CPA initial Sw");
}
}

int main()
{
    try
    {
        {
            std::vector<MPMC::WellPerforation<int>> perforations{{0, 1.0e-12}};
            const auto well = MPMC::cases::makeWell<Indices, int>(
                WellConfig::wells.front(), 0, std::move(perforations));
            require(well.id == WellConfig::wells.front().id,
                    "case well factory must preserve the configured well id");
            require(well.control == MPMC::cases::detail::toNaturalControl(
                        WellConfig::wells.front().control),
                    "case well factory must preserve the configured well control");
        }

        checkPhaseSuppressionCodec();
        checkConfiguration();
        checkInitialFlash();
        checkPrNearPhaseBoundaryHysteresis();
        checkSwNearPhaseBoundaryHysteresis();
        checkConsistentPhaseBoundaryHysteresis();
        checkLatePhaseBoundaryPlateausAcrossEos();
        checkV58CpaMinimumDtPlateauRegression();
        checkStrongOilInstabilityReappears();
        std::cout << "[PASS] three_eos_3d_case_test\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[FAIL] three_eos_3d_case_test: " << error.what() << '\n';
        return 1;
    }
}
