/**
 * @file h2o_co2_nc10_2d_benchmark.cpp
 * @brief 新热力学 H2O-CO2-nC10 规则二维对比算例主程序。
 */
#include <case/petsc_custom_hooks.hpp>
#include "case_config.hpp"
#include "case_fluid.hpp"
#include "well_config.hpp"

#include <case/case_support.hpp>
#include <case/well_factory.hpp>
#include <indices/indices.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
#include <natural/thermo/three_phase_flash.hpp>
#include <structuredgrid/grid_report.hpp>
#include <structuredgrid/structuredgrid.hpp>
#include <well/peaceman.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace Case
{
using Config = CaseConfig::Config;
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents, CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells, CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption, CaseConfig::Model::enableLandTrapping,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ADIndices<ModelConfig>;
using Grid = MPMC::StructuredGridCore;
using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;
using Composition = std::array<double, CaseConfig::CommonFluid::N>;

static_assert(Indices::fullyCompositionalThreePhase);
static_assert(CaseConfig::Grid::nx * CaseConfig::Grid::ny * CaseConfig::Grid::nz == 1200);

enum class EosChoice { Pr, Sw, Cpa };

const char *eosToken(EosChoice choice)
{
    switch (choice) {
    case EosChoice::Pr: return "new-pr";
    case EosChoice::Sw: return "new-sw";
    case EosChoice::Cpa: return "new-cpa";
    }
    return "unknown";
}

EosChoice readEosChoice()
{
    char value[32] = "pr";
    PetscCallAbort(PETSC_COMM_WORLD,
        PetscOptionsGetString(nullptr, nullptr, "-eos", value, sizeof(value), nullptr));
    const std::string token(value);
    if (token == "pr" || token == "PR" || token == "new-pr") return EosChoice::Pr;
    if (token == "sw" || token == "SW" || token == "new-sw") return EosChoice::Sw;
    if (token == "cpa" || token == "CPA" || token == "new-cpa") return EosChoice::Cpa;
    throw std::invalid_argument("-eos must be pr, sw, or cpa.");
}

MPMC::FluidSystem<Indices> makeFluid(EosChoice choice)
{
    MPMC::FluidSystem<Indices> fluid = [&] {
        switch (choice) {
        case EosChoice::Pr:
            return MPMC::cases::makeFluidSystem<Indices, CaseConfig::PrFactoryConfig>();
        case EosChoice::Sw:
            return MPMC::cases::makeFluidSystem<Indices, CaseConfig::SwFactoryConfig>();
        case EosChoice::Cpa:
            return MPMC::cases::makeFluidSystem<Indices, CaseConfig::CpaFactoryConfig>();
        }
        throw std::logic_error("Unknown EOS choice.");
    }();
    BenchmarkCaseFluid::apply(fluid);
    return fluid;
}

struct InitialEquilibrium
{
    Composition overall{};
    MPMC::CubicThreePhaseFlash<Indices>::Result flash{};
};

InitialEquilibrium findInitialEquilibrium(const MPMC::FluidSystem<Indices> &fluid)
{
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = CaseConfig::CommonFluid::waterComponent;
    const MPMC::CubicThreePhaseFlash<Indices> flash(fluid.eos, options);
#ifdef MPMC_H2O_CO2_NC10_FIXED_INITIAL_H2O_MOLE_FRACTION
    // Above the water critical point the fully-compositional EOS models need
    // not retain the low-temperature O+W tie line.  A case may therefore ask
    // for a common, stable overall H2O inventory instead of forcing Sw=0.20.
#ifndef MPMC_H2O_CO2_NC10_FIXED_INITIAL_CO2_MOLE_FRACTION
#define MPMC_H2O_CO2_NC10_FIXED_INITIAL_CO2_MOLE_FRACTION 0.0
#endif
    constexpr double fixedWater =
        MPMC_H2O_CO2_NC10_FIXED_INITIAL_H2O_MOLE_FRACTION;
    constexpr double fixedCo2 =
        MPMC_H2O_CO2_NC10_FIXED_INITIAL_CO2_MOLE_FRACTION;
    static_assert(fixedWater >= 0.0 && fixedCo2 >= 0.0 &&
                  fixedWater + fixedCo2 <= 1.0,
                  "Fixed initial overall composition must be normalized.");
    const Composition fixed{
        fixedWater, fixedCo2, 1.0 - fixedWater - fixedCo2};
    const auto result = flash.flash(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature, fixed);
    if (!result.converged)
        throw std::runtime_error("Fixed initial-composition flash did not converge.");
    return InitialEquilibrium{fixed, result};
#else
    auto evaluate = [&](double zw) {
        const Composition z{zw, 0.0, 1.0 - zw};
        const auto result = flash.flash(
            CaseConfig::InitialState::pressure,
            CaseConfig::InitialState::temperature, z);
        if (!result.converged)
            throw std::runtime_error("H2O-nC10 initial-equilibrium flash did not converge.");
        return InitialEquilibrium{z, result};
    };

    // The target Sw=0.20 corresponds to z_H2O near 0.67 because the aqueous
    // molar density is much larger than the nC10-rich molar density.  Bracket
    // only the physical two-phase neighbourhood; pure-endpoint stability
    // probes are deliberately outside this initialization problem.
    double lower = 0.50;
    double upper = 0.85;
    auto low = evaluate(lower);
    auto high = evaluate(upper);
    const double target = BenchmarkCommon::targetWaterSaturation;
    if (!(low.flash.saturation[Indices::Phase::water] <= target &&
          high.flash.saturation[Indices::Phase::water] >= target))
        throw std::runtime_error("Cannot bracket the target initial water saturation.");

    InitialEquilibrium middle{};
    for (int iteration = 0; iteration < 70; ++iteration)
    {
        const double zw = 0.5 * (lower + upper);
        middle = evaluate(zw);
        if (middle.flash.saturation[Indices::Phase::water] < target)
            lower = zw;
        else
            upper = zw;
    }
    if (std::abs(middle.flash.saturation[Indices::Phase::water] - target) > 1.0e-10)
        throw std::runtime_error("Initial-equilibrium saturation match failed.");
    return middle;
#endif
}

void initializeRock(Grid &grid)
{
    grid.initializeRockProperties();
    auto permeability = grid.vecGetArray<3>(grid.permeabilityVector);
    auto porosity = grid.vecGetArray<1>(grid.porosityVector);
    const auto owned = grid.ownedRegion();
    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                permeability[k][j][i] = {
                    CaseConfig::Rock::kx(i,j,k), CaseConfig::Rock::ky(i,j,k),
                    CaseConfig::Rock::kz(i,j,k)};
                porosity[k][j][i][0] = CaseConfig::Rock::porosity(i,j,k);
            }
    grid.vecRestoreArray<3>(grid.permeabilityVector, permeability);
    grid.vecRestoreArray<1>(grid.porosityVector, porosity);
}

std::vector<Well> makeWells(Grid &grid)
{
    auto wells = MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid, WellConfig::wells,
        [](Grid &g, const auto &definition, int i, int j, int k) {
            const auto size = g.cellSize({i,j,k});
            return MPMC::verticalPeacemanWellIndex({
                size[0], size[1], size[2], CaseConfig::Rock::kx(i,j,k),
                CaseConfig::Rock::ky(i,j,k), definition.radius, definition.skin});
        });
    wells[0].control = MPMC::WellControl::ReservoirTotalRate;
    wells[0].primaryControl = MPMC::WellControl::ReservoirTotalRate;
    return wells;
}

void initializeWellPressures(Grid &grid, Vec solution)
{
    auto values = grid.vecGetArray<Indices::numPrimaryVariables>(solution);
    const auto owned = grid.ownedRegion();
    for (const auto &definition : WellConfig::wells)
    {
        const int i = definition.completion.i;
        const int j = definition.completion.j;
        const int k = definition.completion.kBegin;
        if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
            j >= owned.yStart && j < owned.yStart + owned.yCount &&
            k >= owned.zStart && k < owned.zStart + owned.zCount)
            values[k][j][i][Indices::Primary::wellPressure] = definition.initialBhp;
    }
    grid.vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

void configureBenchmarkSolverDefaults()
{
    PetscBool lineSearchSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(
        nullptr, nullptr, "-snes_linesearch_type", &lineSearchSet));
    if (!lineSearchSet)
        PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(
            nullptr, "-snes_linesearch_type", "basic"));

    PetscBool stepToleranceSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(
        nullptr, nullptr, "-snes_stol", &stepToleranceSet));
    if (!stepToleranceSet)
        PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(
            nullptr, "-snes_stol", "1e-12"));
}

int run()
{
    PetscMPIInt rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPMC::cases::validateCaseConfig<Indices, Config>();
    const auto choice = readEosChoice();
    auto run = MPMC::cases::readRunOptions<Config>();
    char requested[PETSC_MAX_PATH_LEN]{};
    PetscBool resultSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetString(
        nullptr, nullptr, "-result_dir", requested, sizeof(requested), &resultSet));
    if (!resultSet) run.resultDirectory = std::string("./results/") + eosToken(choice);
    MPMC::cases::printRunSummary(CaseConfig::name, run, rank);

    Grid grid(CaseConfig::Grid::nx, CaseConfig::Grid::ny, CaseConfig::Grid::nz,
              MPMC::GridExtent{CaseConfig::Grid::lx, CaseConfig::Grid::ly,
                               CaseConfig::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);

    auto fluid = makeFluid(choice);
    const auto initial = findInitialEquilibrium(fluid);
    if (rank == 0)
        PetscPrintf(PETSC_COMM_SELF,
            "[BENCHMARK] model=%s grid=60x20x1 PV=%.12g m3 rate=%.12g m3/s "
            "inj_max_bhp=%.12g Pa prod_bhp=%.12g Pa\n"
            "[INIT] z(H2O/CO2/nC10)=%.12g/%.12g/%.12g "
            "S(O/G/W)=%.12g/%.12g/%.12g\n",
            eosToken(choice), BenchmarkCommon::poreVolume, BenchmarkCommon::injectionRate,
            BenchmarkCommon::injectorMaximumBhp, BenchmarkCommon::producerBhp,
            initial.overall[0], initial.overall[1], initial.overall[2],
            initial.flash.saturation[0], initial.flash.saturation[1], initial.flash.saturation[2]);

    auto runtimeOptions = MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(run);
    runtimeOptions.scaling.enabled = true;
    // The generic 1 m3/s floor is four orders above this benchmark's rate and
    // would make an appreciable rate error look converged in the scaled norm.
    runtimeOptions.scaling.rateWellResidualFloor = BenchmarkCommon::injectionRate;
    runtimeOptions.maximumAqueousUnsupportedNewtonMoleFraction =
        CaseConfig::CommonFluid::aqueousVolumeMaximumUnsupportedMoleFraction;
    runtimeOptions.aqueousNewtonWaterComponent = CaseConfig::CommonFluid::waterComponent;
    runtimeOptions.aqueousNewtonCo2Component = CaseConfig::CommonFluid::co2Component;
    runtimeOptions.scaling.validate();
    Runtime runtime(grid, fluid, runtimeOptions);
    runtime.setWells(makeWells(grid));
    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    runtime.initializeUniformFromPTZ(solution, CaseConfig::InitialState::pressure,
                                     CaseConfig::InitialState::temperature, initial.overall);
    initializeWellPressures(grid, solution);
    runtime.initializeHistory(solution);
    configureBenchmarkSolverDefaults();
    MPMC::cases::NaturalSolver<Runtime> solver(runtime, grid.dm(Indices::numPrimaryVariables));
    MPMC::cases::runTimeLoop<Indices, Runtime, Config>(
        runtime, solver.snes(), solution, run, rank);
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&solution));
    return 0;
}
} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    PetscCallAbort(PETSC_COMM_WORLD, PetscInitialize(&argc, &argv, nullptr, nullptr));
    int code = 0;
    try { code = Case::run(); }
    catch (const std::exception &error) {
        PetscPrintf(PETSC_COMM_WORLD, "[ERROR][CASE] %s\n", error.what());
        code = 1;
    }
    PetscCallAbort(PETSC_COMM_WORLD, PetscFinalize());
    return code;
}
