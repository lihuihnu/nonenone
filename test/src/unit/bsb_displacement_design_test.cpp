/**
 * @file bsb_displacement_design_test.cpp
 * @brief 对 CO2/SCW–BSB 实验矩阵执行无 PETSc 的配置契约与闪蒸预检。
 */
#include "../../../case/h2o_co2_bsb_lumped_3d_lab/case_config.hpp"
#include "../../../case/h2o_co2_bsb_lumped_3d_lab/well_config.hpp"

#include <case/case_validation.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/state/cell_property_evaluator.hpp>
#include <natural/state/state_codec.hpp>
#include <natural/state/three_phase_equilibrium.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
using ModelConfig = MPMC::CompositionalModelConfig<
    LabScw3D::Model::numberOfComponents,
    LabScw3D::Model::hasWater,
    LabScw3D::Model::hasWells,
    LabScw3D::Model::enableDissolution,
    LabScw3D::Model::enableAdsorption,
    LabScw3D::Model::enableLandTrapping,
    LabScw3D::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance,
          const std::string &message)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (std::abs(actual - expected) > tolerance * scale)
        throw std::runtime_error(message);
}

MPMC::CubicThreePhaseFlash<Indices>::Result flash(
    const MPMC::CubicThreePhaseFlash<Indices> &solver,
    const std::array<double, Indices::numComponents> &composition)
{
    return solver.flash(
        LabScw3D::InitialState::pressure,
        LabScw3D::InitialState::temperature,
        composition);
}

double maximumClosureError(
    const MPMC::CubicThreePhaseFlash<Indices>::Result &result,
    const std::array<double, Indices::numComponents> &feed)
{
    double maximum = 0.0;
    for (int component = 0; component < Indices::numComponents; ++component)
    {
        double reconstructed = 0.0;
        for (int phase = 0; phase < Indices::numPhases; ++phase)
            reconstructed += result.phaseMoleFraction[static_cast<std::size_t>(phase)] *
                result.composition[static_cast<std::size_t>(phase)]
                                  [static_cast<std::size_t>(component)];
        maximum = std::max(maximum,
            std::abs(reconstructed - feed[static_cast<std::size_t>(component)]));
    }
    return maximum;
}

void checkConfiguration()
{
    static_assert(Indices::fullyCompositionalThreePhase);
    static_assert(LabScw3D::Grid::nx == 60);
    static_assert(LabScw3D::Grid::ny == 20);
    static_assert(LabScw3D::Grid::nz == 1);
    static_assert(LabScw3D::Grid::lx == 300.0);
    static_assert(LabScw3D::Grid::ly == 100.0);
    static_assert(LabScw3D::Grid::lz == 5.0);
    static_assert(LabScw3D::CommonFluid::temperature == 653.15);
    static_assert(LabScw3D::InitialState::pressure == 28.0e6);
    static_assert(LabScw3D::CommonFluid::useIapwsGarciaAqueousVolume);
    static_assert(LabScw3D::CommonFluid::useMcBrideWrightAqueousViscosity);

    MPMC::cases::validateCaseConfig<Indices, LabScw3D::Config>();
    near(BenchmarkCommon::poreVolume, 30000.0, 1.0e-14,
         "shared nC10 pore volume");
    near(BenchmarkCommon::injectionRate * BenchmarkCommon::year /
             BenchmarkCommon::poreVolume,
         0.10, 1.0e-14, "shared nC10 injection PVI per year");

    double compositionSum = 0.0;
    for (double value : LabScw3D::InitialState::overallComposition)
        compositionSum += value;
    near(compositionSum, 1.0, 1.0e-14, "initial composition normalization");

    const auto wells = LabScw3DWell::definitions(
        60, 20, 1, BenchmarkCommon::injectionRate);
    require(wells[0].completion.i == 0 && wells[0].completion.j == 9,
            "injector must match the nC10 grid");
    require(wells[1].completion.i == 59 && wells[1].completion.j == 9,
            "producer must match the nC10 grid");
    near(wells[1].target, BenchmarkCommon::producerBhp, 0.0,
         "producer BHP must match the nC10 SCW case");
}

void checkPrFlashPreflight()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, LabScw3D::PrFactoryConfig>();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = LabScw3D::CommonFluid::waterComponent;
    const MPMC::CubicThreePhaseFlash<Indices> solver(fluid.eos, options);

    const auto initial = flash(solver, LabScw3D::InitialState::overallComposition);
    require(initial.converged, "initial BSB P-T-z flash must converge");
    require(initial.phaseMoleFraction[0] > 0.0,
            "initial BSB state must contain an oil-rich phase");
    require(initial.presence.count() == 1 && initial.phaseMoleFraction[0] == 1.0,
            "initial BSB preflight is expected to be an oil-only state");
    require(maximumClosureError(initial, LabScw3D::InitialState::overallComposition) <= 1.0e-10,
            "initial flash composition closure");

    MPMC::FullyCompositionalThreePhaseEquilibrium<Indices> equilibrium(fluid);
    std::array<double, Indices::numPrimaryVariables> primary{};
    primary[Indices::Primary::pressure] = LabScw3D::InitialState::pressure;
    MPMC::PhaseStateData<Indices> phaseState;
    equilibrium.assignFlashResult(primary, phaseState, initial);
    const auto state = MPMC::CellStateCodec<Indices>::decode(primary, phaseState);
    MPMC::CellPropertyEvaluator<Indices> evaluator(fluid);
    const auto properties = evaluator.evaluate(state, BenchmarkCommon::Rock::porosityValue);
    require(std::isfinite(properties.density[Indices::Phase::water]) &&
                properties.density[Indices::Phase::water] > 0.0,
            "inactive water reference density must remain finite");
    near(properties.mobility[Indices::Phase::water], 0.0, 0.0,
         "inactive water phase mobility");

    constexpr std::array<double, 5> scwLevels{0.0, 0.25, 0.50, 0.75, 1.0};
    for (double scw : scwLevels)
    {
        auto mixed = LabScw3D::InitialState::overallComposition;
        for (double &value : mixed)
            value *= 0.99;
        mixed[0] += 0.01 * scw;
        mixed[1] += 0.01 * (1.0 - scw);
        const auto result = flash(solver, mixed);
        require(result.converged, "1% inlet-mixed flash must converge");
        require(maximumClosureError(result, mixed) <= 1.0e-10,
                "1% inlet-mixed flash composition closure");
    }
}

void checkHydrocarbonRichWaterLabelRegression()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, LabScw3D::PrFactoryConfig>();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = LabScw3D::CommonFluid::waterComponent;
    const MPMC::CubicThreePhaseFlash<Indices> solver(fluid.eos, options);
    const std::array<double, Indices::numComponents> composition{
        0.57589755842924628, 0.17086252972424726,
        0.10574594269902246, 0.08658849930051303,
        0.042219755560873616, 0.018685714286097353};
    constexpr double pressure = 2.8119450976194970e7;
    MPMC::PhasePresence oilOnly = MPMC::PhasePresence::oilOnly();
    const std::array<std::array<double, Indices::numComponents>, 3> phases{
        composition, composition, composition};
    const auto stability = solver.stabilityTest(
        pressure,
        LabScw3D::InitialState::temperature,
        composition,
        oilOnly,
        phases);
    require(stability.valid, "phase-boundary stability test must be valid");
    require(stability.missingPhaseUnstable[1] &&
            !stability.missingPhaseUnstable[2],
            "only the admissible PR gas trial may enter the active set");
    require(!fluid.eos.aqueousVolumeCompositionSupported(
                stability.incipientComposition[2]),
            "incipient hydrocarbon-rich phase must remain outside aqueous closure");

    const auto result = solver.flash(
        pressure,
        LabScw3D::InitialState::temperature,
        composition);
    require(result.converged, "phase-boundary unrestricted flash must converge");
    require(result.presence.contains(MPMC::CompositionalPhase::Oil) &&
            result.presence.contains(MPMC::CompositionalPhase::Gas) &&
            !result.presence.contains(MPMC::CompositionalPhase::Water),
            "hydrocarbon-rich split must remain an O+G state, not an aqueous state");
    require(maximumClosureError(result, composition) <= 1.0e-10,
            "phase-boundary flash composition closure");
}

} // namespace

int main()
{
    try
    {
        checkConfiguration();
        checkPrFlashPreflight();
        checkHydrocarbonRichWaterLabelRegression();
        std::cout << "[PASS] bsb_displacement_design_test\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[FAIL] bsb_displacement_design_test: " << error.what() << '\n';
        return 1;
    }
}
