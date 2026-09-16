/**
 * @file h2o_co2_bsb_lumped_3d_lab.cpp
 * @brief 在共用 nC10 二维网格上运行五点 CO2/SCW 驱替对照。
 */
#include <case/petsc_custom_hooks.hpp>
#include "case_config.hpp"
#include "well_config.hpp"

#include <case/case_support.hpp>
#include <case/natural_scaling_options.hpp>
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

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace Case
{

using Config = LabScw3D::Config;
using ModelConfig = MPMC::CompositionalModelConfig<
    Config::Model::numberOfComponents, Config::Model::hasWater,
    Config::Model::hasWells, Config::Model::enableDissolution,
    Config::Model::enableAdsorption, Config::Model::enableLandTrapping,
    Config::Model::phaseBehavior>;
using Indices = MPMC::ADIndices<ModelConfig>;
using Grid = MPMC::StructuredGridCore;
using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;

static_assert(Indices::fullyCompositionalThreePhase);
static_assert(LabScw3D::CommonFluid::temperature > 647.096);
static_assert(LabScw3D::InitialState::pressure > 22.064e6);

enum class EosChoice { Pr, Sw };

struct StudyOptions
{
    EosChoice eos{EosChoice::Pr};
    double xH2OFeed{0.50};
    int nx{LabScw3D::Grid::nx};
    int ny{LabScw3D::Grid::ny};
    int nz{LabScw3D::Grid::nz};
    double targetPvi{0.10};
    double dtPvi{0.002};
};

const char *eosToken(EosChoice eos)
{
    return eos == EosChoice::Pr ? "pr" : "sw";
}

StudyOptions readStudyOptions()
{
    char eosText[16] = "";
    PetscBool eosSet = PETSC_FALSE;
    PetscBool feedSet = PETSC_FALSE;
    PetscInt nx = LabScw3D::Grid::nx;
    PetscInt ny = LabScw3D::Grid::ny;
    PetscInt nz = LabScw3D::Grid::nz;
    PetscReal xH2O = 0.50;
    PetscReal targetPvi = 0.10;
    PetscReal dtPvi = 0.002;

    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetString(
        nullptr, nullptr, "-eos", eosText, sizeof(eosText), &eosSet));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetReal(
        nullptr, nullptr, "-x_h2o_feed", &xH2O, &feedSet));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetInt(nullptr, nullptr, "-nx", &nx, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetInt(nullptr, nullptr, "-ny", &ny, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetInt(nullptr, nullptr, "-nz", &nz, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetReal(
        nullptr, nullptr, "-target_pv", &targetPvi, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetReal(
        nullptr, nullptr, "-dt_pv", &dtPvi, nullptr));

    if (!eosSet)
        throw std::invalid_argument("-eos pr|sw is required.");
    if (!feedSet)
        throw std::invalid_argument("-x_h2o_feed in [0,1] is required.");

    StudyOptions result;
    const std::string eos(eosText);
    if (eos == "pr" || eos == "PR")
        result.eos = EosChoice::Pr;
    else if (eos == "sw" || eos == "SW")
        result.eos = EosChoice::Sw;
    else
        throw std::invalid_argument("-eos must be pr or sw.");
    result.xH2OFeed = static_cast<double>(xH2O);
    result.nx = static_cast<int>(nx);
    result.ny = static_cast<int>(ny);
    result.nz = static_cast<int>(nz);
    result.targetPvi = static_cast<double>(targetPvi);
    result.dtPvi = static_cast<double>(dtPvi);

    if (!std::isfinite(result.xH2OFeed) || result.xH2OFeed < 0.0 || result.xH2OFeed > 1.0)
        throw std::invalid_argument("-x_h2o_feed must be finite and in [0,1].");
    if (result.nx != LabScw3D::Grid::nx ||
        result.ny != LabScw3D::Grid::ny ||
        result.nz != LabScw3D::Grid::nz)
        throw std::invalid_argument(
            "The primary comparison fixes the H2O-CO2-nC10 grid at 60x20x1.");
    if (!(result.targetPvi > 0.0) || !(result.dtPvi > 0.0) ||
        result.dtPvi > result.targetPvi)
        throw std::invalid_argument("PVI horizon and step must be positive, with dt <= target.");
    return result;
}

MPMC::FluidSystem<Indices> makeFluid(EosChoice eos)
{
    auto fluid = eos == EosChoice::Pr
        ? MPMC::cases::makeFluidSystem<Indices, LabScw3D::PrFactoryConfig>()
        : MPMC::cases::makeFluidSystem<Indices, LabScw3D::SwFactoryConfig>();
    BenchmarkCommon::applyCommonRelativePermeability(fluid);
    return fluid;
}

void initializeRock(Grid &grid)
{
    grid.initializeRockProperties();
    auto rock = grid.rockWriteView();
    const auto owned = grid.ownedRegion();
    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                rock.permeability(i, j, k) = {
                    LabScw3D::Rock::kx(i, j, k), LabScw3D::Rock::ky(i, j, k),
                    LabScw3D::Rock::kz(i, j, k)};
                rock.porosity(i, j, k) = LabScw3D::Rock::porosity(i, j, k);
            }
}

template <class Definition>
double wellIndex(Grid &grid, const Definition &definition, int i, int j, int k)
{
    const auto size = grid.cellSize({i, j, k});
    return MPMC::verticalPeacemanWellIndex({
        size[0], size[1], size[2], LabScw3D::Rock::kx(i, j, k),
        LabScw3D::Rock::ky(i, j, k), definition.radius, definition.skin});
}

struct InjectionState
{
    std::array<double, Indices::numPhases> phaseSlotFraction{};
    std::array<double, Indices::numComponents> componentMassFraction{};
};

InjectionState makeInjectionState(double xH2O)
{
    InjectionState state;
    const double mh = LabScw3D::CommonFluid::molarMass[0];
    const double mc = LabScw3D::CommonFluid::molarMass[1];
    const double feedMass = xH2O * mh + (1.0 - xH2O) * mc;
    state.componentMassFraction[0] = xH2O * mh / feedMass;
    state.componentMassFraction[1] = (1.0 - xH2O) * mc / feedMass;

    state.phaseSlotFraction[Indices::Phase::water] = xH2O;
    state.phaseSlotFraction[Indices::Phase::vapor] = 1.0 - xH2O;
    return state;
}

std::vector<Well> makeWells(
    Grid &grid,
    const std::array<LabScw3DWell::WellDefinition, 2> &definitions,
    const InjectionState &injection)
{
    auto wells = MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid, definitions,
        [](Grid &g, const auto &definition, int i, int j, int k) {
            return wellIndex(g, definition, i, j, k);
        });
    auto &injector = wells[0];
    injector.control = MPMC::WellControl::ReservoirTotalRate;
    injector.primaryControl = MPMC::WellControl::ReservoirTotalRate;
    injector.injectionPhaseFraction = injection.phaseSlotFraction;
    injector.injectionComponentMassFraction = injection.componentMassFraction;
    return wells;
}

void initializeWellPressureGuesses(
    Grid &grid, Vec solution,
    const std::array<LabScw3DWell::WellDefinition, 2> &definitions)
{
    auto values = grid.template vecGetArray<Indices::numPrimaryVariables>(solution);
    const auto owned = grid.ownedRegion();
    for (const auto &definition : definitions)
    {
        const int i = definition.completion.i;
        const int j = definition.completion.j;
        const int k = definition.completion.kBegin;
        if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
            j >= owned.yStart && j < owned.yStart + owned.yCount &&
            k >= owned.zStart && k < owned.zStart + owned.zCount)
            values[k][j][i][Indices::Primary::wellPressure] = definition.initialBhp;
    }
    grid.template vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

double initialReferenceDensity(const MPMC::FluidSystem<Indices> &fluid)
{
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = LabScw3D::CommonFluid::waterComponent;
    const MPMC::CubicThreePhaseFlash<Indices> flash(fluid.eos, options);
    const auto result = flash.flash(
        LabScw3D::InitialState::pressure,
        LabScw3D::InitialState::temperature,
        LabScw3D::InitialState::overallComposition);
    if (!result.converged)
        throw std::runtime_error("Initial P-T-z flash did not converge.");

    const auto phase = static_cast<std::size_t>(
        std::distance(result.phaseMoleFraction.begin(),
                      std::max_element(result.phaseMoleFraction.begin(),
                                       result.phaseMoleFraction.end())));
    double meanMolarMass = 0.0;
    for (std::size_t c = 0; c < LabScw3D::CommonFluid::N; ++c)
        meanMolarMass += result.composition[phase][c] *
                         LabScw3D::CommonFluid::molarMass[c];
    const double density = result.molarDensity[phase] * meanMolarMass;
    if (!(density > 0.0) || !std::isfinite(density))
        throw std::runtime_error("Initial flash returned an invalid reference density.");
    return density;
}

void applyHydrostaticPressure(Grid &grid, Vec solution, double density)
{
    auto values = grid.template vecGetArray<Indices::numPrimaryVariables>(solution);
    const auto owned = grid.ownedRegion();
    const auto dims = grid.dimensions();
    const double dz = LabScw3D::Grid::lz / static_cast<double>(dims[2]);
    const double topCellCenter = LabScw3D::Grid::lz - 0.5 * dz;
    constexpr double gravity = 9.80665;
    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
    {
        const double z = (static_cast<double>(k) + 0.5) * dz;
        const double pressure = LabScw3D::InitialState::pressure +
                                density * gravity * (topCellCenter - z);
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
                values[k][j][i][Indices::Primary::pressure] = pressure;
    }
    grid.template vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

void writeRunMetadata(
    const MPMC::cases::RunOptions &run, const StudyOptions &study,
    const InjectionState &injection, double effectiveDtPvi, int rank)
{
    if (rank != 0)
        return;
    std::filesystem::create_directories(run.resultDirectory);
    std::ofstream stream(std::filesystem::path(run.resultDirectory) / "design_metadata.csv");
    if (!stream)
        throw std::runtime_error("Failed to write design_metadata.csv.");
    stream << std::setprecision(17)
           << "eos,x_scw_feed,x_co2_feed,y_h2o_feed,y_co2_feed,nx,ny,nz,"
              "target_pvi,dt_pvi,total_reservoir_rate_m3_s,pore_volume_m3,"
              "temperature_K,initial_pressure_Pa,producer_bhp_Pa\n"
           << eosToken(study.eos) << ',' << study.xH2OFeed << ','
           << (1.0 - study.xH2OFeed) << ','
           << injection.componentMassFraction[0] << ','
           << injection.componentMassFraction[1] << ','
           << study.nx << ',' << study.ny << ',' << study.nz << ','
           << study.targetPvi << ',' << effectiveDtPvi << ','
           << BenchmarkCommon::injectionRate << ',' << BenchmarkCommon::poreVolume << ','
           << LabScw3D::CommonFluid::temperature << ','
           << LabScw3D::InitialState::pressure << ',' << BenchmarkCommon::producerBhp << '\n';
}

void configureSolverDefaults()
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
    const StudyOptions study = readStudyOptions();
    constexpr double totalRate = BenchmarkCommon::injectionRate;
    constexpr double poreVolume = BenchmarkCommon::poreVolume;

    auto runOptions = MPMC::cases::readRunOptions<Config>();
    const std::size_t intervals = static_cast<std::size_t>(
        std::ceil(study.targetPvi / study.dtPvi));
    const double effectiveDtPvi = study.targetPvi / static_cast<double>(intervals);
    runOptions.numberOfSteps = intervals;
    runOptions.dtDays = (poreVolume / totalRate) * effectiveDtPvi /
                        LabScw3D::secondsPerDay;

    PetscBool resultSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(
        nullptr, nullptr, "-result_dir", &resultSet));
    if (!resultSet)
    {
        const int percent = static_cast<int>(std::lround(100.0 * study.xH2OFeed));
        const int co2Percent = 100 - percent;
        runOptions.resultDirectory = std::string("./results/") + eosToken(study.eos) +
            "_CO2_" + std::to_string(co2Percent) + "_SCW_" +
            std::to_string(percent);
    }
    MPMC::cases::printRunSummary(Config::name, runOptions, rank);

    Grid grid(study.nx, study.ny, study.nz,
              MPMC::GridExtent{LabScw3D::Grid::lx, LabScw3D::Grid::ly,
                               LabScw3D::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);

    auto fluid = makeFluid(study.eos);
    const auto injection = makeInjectionState(study.xH2OFeed);
    const auto definitions = LabScw3DWell::definitions(
        study.nx, study.ny, study.nz, totalRate);
    writeRunMetadata(runOptions, study, injection, effectiveDtPvi, rank);

    auto runtimeOptions =
        MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(runOptions);
    MPMC::cases::applyNaturalScalingPetscOptions(runtimeOptions, rank);
    runtimeOptions.scaling.enabled = true;
    runtimeOptions.scaling.rateWellResidualFloor = totalRate;
    runtimeOptions.scaling.validate();

    Runtime runtime(grid, fluid, runtimeOptions);
    runtime.setWells(makeWells(grid, definitions, injection));
    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    runtime.initializeUniformFromPTZ(
        solution, LabScw3D::InitialState::pressure,
        LabScw3D::InitialState::temperature,
        LabScw3D::InitialState::overallComposition);
    applyHydrostaticPressure(grid, solution, initialReferenceDensity(fluid));
    initializeWellPressureGuesses(grid, solution, definitions);
    runtime.updateState(solution);
    runtime.initializeHistory(solution);

    configureSolverDefaults();
    MPMC::cases::NaturalSolver<Runtime> solver(
        runtime, grid.dm(Indices::numPrimaryVariables));
    MPMC::cases::runTimeLoop<Indices, Runtime, Config>(
        runtime, solver.snes(), solution, runOptions, rank);
    PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&solution));
    return 0;
}

} // namespace Case

MPMC_DEFINE_NATURAL_PETSC_CUSTOM_HOOKS(Case::Runtime)

int main(int argc, char **argv)
{
    PetscCallAbort(PETSC_COMM_WORLD, PetscInitialize(&argc, &argv, nullptr, nullptr));
    int code = 0;
    try
    {
        code = Case::run();
    }
    catch (const std::exception &error)
    {
        PetscPrintf(PETSC_COMM_WORLD, "[ERROR][CASE] %s\n", error.what());
        code = 1;
    }
    PetscCallAbort(PETSC_COMM_WORLD, PetscFinalize());
    return code;
}
