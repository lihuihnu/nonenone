/**
 * @file sun2024_scw_co2_nc16_factorial_1d.cpp
 * @brief 可复现Sun-2024超临界水-CO2二乘二实验及等流量对照。
 */
#include <case/petsc_custom_hooks.hpp>
#include "case_config.hpp"
#include "case_fluid.hpp"
#include "experiment_matrix.hpp"
#include "well_config.hpp"

#include <case/case_support.hpp>
#include <case/well_factory.hpp>
#include <indices/indices.hpp>
#include <natural/petsc/natural_structuredgrid_runtime.hpp>
#include <natural/thermo/aqueous_volume.hpp>
#include <structuredgrid/grid_report.hpp>
#include <structuredgrid/structuredgrid.hpp>
#include <well/peaceman.hpp>
#include <well/well.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
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
using ScalarIndices = MPMC::ScalarIndices<ModelConfig>;
using Grid = MPMC::StructuredGridCore;
using Runtime = MPMC::NaturalStructuredGridRuntime<Indices, Grid>;
using Well = MPMC::NaturalWell<Indices>;

static_assert(!Indices::fullyCompositionalThreePhase);
static_assert(Indices::hasIndependentWaterConservation);
static_assert(Sun2024Factorial::temperature >
              Sun2024Factorial::waterCriticalTemperature);

struct StudyOptions
{
    const Sun2024Factorial::ExperimentDefinition *experiment{nullptr};
    int nx{48};
    double targetPoreVolumes{Sun2024Factorial::defaultTargetPoreVolumes};
    double outputStepPoreVolumes{
        Sun2024Factorial::defaultOutputStepPoreVolumes};
    Sun2024Factorial::RateBasis rateBasis{
        Sun2024Factorial::RateBasis::InSituVolume};
    Sun2024Factorial::ReferenceDensities referenceDensities;
};

StudyOptions readStudyOptions()
{
    char experimentId[32]{};
    char rateBasis[64]{};
    PetscBool experimentSet = PETSC_FALSE;
    PetscBool rateBasisSet = PETSC_FALSE;
    PetscInt nx = 48;
    PetscReal targetPv = Sun2024Factorial::defaultTargetPoreVolumes;
    PetscReal dtPv = Sun2024Factorial::defaultOutputStepPoreVolumes;
    PetscReal waterReferenceDensity =
        std::numeric_limits<PetscReal>::quiet_NaN();
    PetscReal co2ReferenceDensity =
        std::numeric_limits<PetscReal>::quiet_NaN();
    PetscBool waterReferenceDensitySet = PETSC_FALSE;
    PetscBool co2ReferenceDensitySet = PETSC_FALSE;

    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetString(
        nullptr, nullptr, "-experiment_id", experimentId,
        sizeof(experimentId), &experimentSet));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetString(
        nullptr, nullptr, "-rate_basis", rateBasis,
        sizeof(rateBasis), &rateBasisSet));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetInt(
        nullptr, nullptr, "-nx", &nx, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetReal(
        nullptr, nullptr, "-target_pv", &targetPv, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetReal(
        nullptr, nullptr, "-dt_pv", &dtPv, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetReal(
        nullptr, nullptr, "-water_reference_density_kg_m3",
        &waterReferenceDensity, &waterReferenceDensitySet));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetReal(
        nullptr, nullptr, "-co2_reference_density_kg_m3",
        &co2ReferenceDensity, &co2ReferenceDensitySet));

    if (!experimentSet)
        throw std::invalid_argument(
            "-experiment_id is required (R04/R06/R10/R12/C23/C24).");
    if (!rateBasisSet)
        throw std::invalid_argument(
            "-rate_basis is required. Use in_situ_volume for an explicitly "
            "declared sensitivity case, or reference_density with measured "
            "phase densities from the experimental flow-meter state.");
    if (nx < 2)
        throw std::invalid_argument("-nx must be at least 2.");
    if (!(targetPv > 0.0) || !std::isfinite(targetPv))
        throw std::invalid_argument("-target_pv must be finite and positive.");
    if (!(dtPv > 0.0) || !std::isfinite(dtPv) || dtPv > targetPv)
        throw std::invalid_argument(
            "-dt_pv must be finite, positive, and no larger than -target_pv.");

    PetscBool genericDtSet = PETSC_FALSE;
    PetscBool genericStepsSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(
        nullptr, nullptr, "-dt", &genericDtSet));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(
        nullptr, nullptr, "-numSteps", &genericStepsSet));
    if (genericDtSet || genericStepsSet)
        throw std::invalid_argument(
            "This case is PVI-controlled. Use -target_pv and -dt_pv, not "
            "-numSteps or -dt.");

    StudyOptions result;
    result.experiment = &Sun2024Factorial::definition(experimentId);
    result.nx = static_cast<int>(nx);
    result.targetPoreVolumes = static_cast<double>(targetPv);
    result.outputStepPoreVolumes = static_cast<double>(dtPv);
    result.rateBasis = Sun2024Factorial::parseRateBasis(rateBasis);
    if (waterReferenceDensitySet)
        result.referenceDensities.waterKgPerM3 = waterReferenceDensity;
    if (co2ReferenceDensitySet)
        result.referenceDensities.co2KgPerM3 = co2ReferenceDensity;
    return result;
}

template <class TargetIndices = Indices>
MPMC::FluidSystem<TargetIndices> makeFluid()
{
    auto fluid = MPMC::cases::makeFluidSystem<
        TargetIndices, CaseConfig::PrFactoryConfig>();
    Sun2024FactorialFluid::apply(fluid);
    return fluid;
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
                    CaseConfig::Rock::kx(i, j, k),
                    CaseConfig::Rock::ky(i, j, k),
                    CaseConfig::Rock::kz(i, j, k)};
                porosity[k][j][i][0] =
                    CaseConfig::Rock::porosity(i, j, k);
            }
    grid.vecRestoreArray<3>(grid.permeabilityVector, permeability);
    grid.vecRestoreArray<1>(grid.porosityVector, porosity);
}

std::array<WellConfig::WellDefinition, 2> makeWellDefinitions(
    const Sun2024Factorial::ExperimentDefinition &experiment,
    const Sun2024Factorial::ResolvedInjection &injection,
    int nx)
{
    auto definitions = WellConfig::wells;
    const double totalRate = injection.totalReservoirRateM3PerS();

    definitions[0].target = totalRate;
    definitions[0].initialBhp = experiment.pressurePa;
    // Figure 1 reports a 35 MPa generator limit; Figure 4 shows multi-MPa
    // inlet-outlet pressure differences.  Guard the apparatus limit instead
    // of clamping the physical pressure transient at 0.2 MPa.
    definitions[0].maximumBhp = 35.0e6;
    definitions[0].completion.i = 0;

    definitions[1].target = experiment.pressurePa;
    definitions[1].initialBhp = experiment.pressurePa;
    definitions[1].completion.i = nx - 1;
    return definitions;
}

std::vector<Well> makeWells(
    Grid &grid,
    const std::array<WellConfig::WellDefinition, 2> &definitions,
    const Sun2024Factorial::ResolvedInjection &injection)
{
    auto wells = MPMC::cases::makeStructuredWells<Indices, PetscInt>(
        grid, definitions,
        [](Grid &g, const auto &definition, int i, int j, int k) {
            const auto size = g.cellSize({i, j, k});
            return MPMC::verticalPeacemanWellIndex({
                size[0], size[1], size[2], CaseConfig::Rock::kx(i, j, k),
                CaseConfig::Rock::ky(i, j, k), definition.radius,
                definition.skin});
        });

    auto &injector = wells[0];
    injector.control = MPMC::WellControl::ReservoirTotalRate;
    injector.primaryControl = MPMC::WellControl::ReservoirTotalRate;
    injector.injectionPhaseFraction.fill(0.0);
    injector.injectionPhaseFraction[Indices::Phase::water] =
        injection.waterReservoirVolumeFraction();
    injector.injectionPhaseFraction[Indices::Phase::vapor] =
        injection.co2ReservoirVolumeFraction();
    injector.injectionComponentMassFraction.fill(0.0);
    injector.injectionComponentMassFraction[
        CaseConfig::Fluid::co2Component] = 1.0;
    return wells;
}

void initializeSolution(
    Grid &grid,
    Vec solution,
    const std::array<WellConfig::WellDefinition, 2> &definitions,
    double pressure)
{
    PetscCallAbort(PETSC_COMM_WORLD, VecSet(solution, 0.0));
    auto values = grid.vecGetArray<Indices::numPrimaryVariables>(solution);
    const auto owned = grid.ownedRegion();
    for (int k = owned.zStart; k < owned.zStart + owned.zCount; ++k)
        for (int j = owned.yStart; j < owned.yStart + owned.yCount; ++j)
            for (int i = owned.xStart; i < owned.xStart + owned.xCount; ++i)
            {
                auto &primary = values[k][j][i];
                primary.fill(0.0);
                primary[Indices::Primary::pressure] = pressure;
                primary[Indices::Primary::liquidComposition[0]] =
                    CaseConfig::InitialState::oilComposition[0];
                primary[Indices::Primary::vaporComposition[0]] =
                    CaseConfig::InitialState::gasComposition[0];
                primary[Indices::Primary::liquidSaturation] =
                    CaseConfig::InitialState::oilSaturation;
                primary[Indices::Primary::vaporSaturation] =
                    CaseConfig::InitialState::gasSaturation;
                primary[Indices::Primary::waterSaturation] =
                    CaseConfig::InitialState::waterSaturation;
            }

    for (const auto &definition : definitions)
    {
        const int i = definition.completion.i;
        const int j = definition.completion.j;
        const int k = definition.completion.kBegin;
        if (i >= owned.xStart && i < owned.xStart + owned.xCount &&
            j >= owned.yStart && j < owned.yStart + owned.yCount &&
            k >= owned.zStart && k < owned.zStart + owned.zCount)
            values[k][j][i][Indices::Primary::wellPressure] =
                definition.initialBhp;
    }
    grid.vecRestoreArray<Indices::numPrimaryVariables>(solution, values);
}

void writeExperimentMetadata(
    const MPMC::cases::RunOptions &run,
    const StudyOptions &study,
    const Sun2024Factorial::ResolvedInjection &injection,
    const Sun2024Factorial::ReservoirDensities &reservoirDensities,
    double effectiveDtPv,
    int rank)
{
    if (rank != 0)
        return;
    std::filesystem::create_directories(run.resultDirectory);
    const auto file = std::filesystem::path(run.resultDirectory) /
        "experiment_metadata.csv";
    std::ofstream stream(file, std::ios::out | std::ios::trunc);
    if (!stream)
        throw std::runtime_error(
            "Failed to create experiment metadata: " + file.string());
    const auto &e = *study.experiment;
    stream << std::setprecision(17)
           << "experiment_id,source_experiment,role,temperature_K,pressure_Pa,"
              "water_rate_reported_mL_min,co2_rate_reported_mL_min,rate_basis,"
              "water_reference_density_kg_m3,co2_reference_density_kg_m3,"
              "water_reservoir_density_kg_m3,co2_reservoir_density_kg_m3,"
              "water_reference_equivalent_mL_min,"
              "co2_reference_equivalent_mL_min,water_mass_rate_kg_s,"
              "co2_mass_rate_kg_s,water_reservoir_rate_m3_s,"
              "co2_reservoir_rate_m3_s,reservoir_total_rate_m3_s,"
              "matched_physical_experiment,pore_volume_m3,target_pv,dt_pv,nx\n"
           << e.id << ',' << e.sourceExperiment << ','
           << Sun2024Factorial::roleName(e.role) << ','
           << Sun2024Factorial::temperature << ',' << e.pressurePa << ','
           << e.waterRateMlPerMinute << ',' << e.co2RateMlPerMinute << ','
           << Sun2024Factorial::rateBasisName(study.rateBasis) << ','
           << study.referenceDensities.waterKgPerM3 << ','
           << study.referenceDensities.co2KgPerM3 << ','
           << reservoirDensities.waterKgPerM3 << ','
           << reservoirDensities.co2KgPerM3 << ','
           << injection.waterReferenceEquivalentMlPerMinute << ','
           << injection.co2ReferenceEquivalentMlPerMinute << ','
           << injection.waterMassRateKgPerS << ','
           << injection.co2MassRateKgPerS << ','
           << injection.waterReservoirRateM3PerS << ','
           << injection.co2ReservoirRateM3PerS << ','
           << injection.totalReservoirRateM3PerS() << ','
           << e.matchedPhysicalExperiment << ','
           << Sun2024Exp12::poreVolume << ',' << study.targetPoreVolumes << ','
           << effectiveDtPv << ',' << study.nx << '\n';
}

void configureSolverDefaults()
{
    PetscBool absoluteToleranceSet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(
        nullptr, nullptr, "-snes_atol", &absoluteToleranceSet));
    if (!absoluteToleranceSet)
        PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsSetValue(
            nullptr, "-snes_atol", "1e-12"));

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
    const auto &experiment = *study.experiment;

    if (!(experiment.pressurePa > Sun2024Factorial::waterCriticalPressure))
        throw std::invalid_argument(
            "Selected pressure is not above the pure-water critical pressure.");

    auto scalarFluid = makeFluid<ScalarIndices>();
    const double waterDensity = MPMC::IapwsIf97WaterDensity::density(
        experiment.pressurePa, Sun2024Factorial::temperature);
    typename MPMC::FluidSystem<ScalarIndices>::Composition pureCo2{};
    pureCo2[CaseConfig::Fluid::co2Component] = 1.0;
    const auto co2Phase = scalarFluid.eos.phaseResult(
        experiment.pressurePa, Sun2024Factorial::temperature, pureCo2, false);
    const auto [co2MassFraction, co2Density, co2Viscosity] =
        scalarFluid.eosFlowProperties(
            experiment.pressurePa, pureCo2, co2Phase.compressibility,
            Sun2024Factorial::temperature, MPMC::CompositionalPhase::Gas);
    (void)co2MassFraction;
    (void)co2Viscosity;
    const Sun2024Factorial::ReservoirDensities reservoirDensities{
        waterDensity, co2Density};
    const auto injection = Sun2024Factorial::resolveInjection(
        experiment, study.rateBasis, study.referenceDensities,
        reservoirDensities);
    auto fluid = makeFluid();

    auto runOptions = MPMC::cases::readRunOptions<Config>();
    const std::size_t intervals = static_cast<std::size_t>(
        std::ceil(study.targetPoreVolumes / study.outputStepPoreVolumes));
    const double effectiveDtPv = study.targetPoreVolumes /
        static_cast<double>(intervals);
    runOptions.numberOfSteps = intervals;
    runOptions.dtDays = Sun2024Factorial::durationDays(
        injection, Sun2024Exp12::poreVolume, effectiveDtPv);

    PetscBool resultDirectorySet = PETSC_FALSE;
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsHasName(
        nullptr, nullptr, "-result_dir", &resultDirectorySet));
    if (!resultDirectorySet)
        runOptions.resultDirectory =
            std::string("./results/factorial/") + experiment.id +
            "_nx" + std::to_string(study.nx);

    MPMC::cases::printRunSummary(CaseConfig::name, runOptions, rank);
    writeExperimentMetadata(
        runOptions, study, injection, reservoirDensities, effectiveDtPv, rank);

    Grid grid(study.nx, 1, 1,
              MPMC::GridExtent{CaseConfig::Grid::lx,
                               CaseConfig::Grid::ly,
                               CaseConfig::Grid::lz});
    grid.setup();
    initializeRock(grid);
    MPMC::printStructuredGridSummary(grid);

    const auto definitions = makeWellDefinitions(
        experiment, injection, study.nx);
    const double totalRate = injection.totalReservoirRateM3PerS();
    if (rank == 0)
        PetscPrintf(PETSC_COMM_SELF,
            "[EXPERIMENT] id=%s source_exp=%d role=%s rate_basis=%s\n"
            "[DOMAIN] grid=%dx1x1 T=%.8g K P=%.8g Pa PV=%.12g m3 "
            "target=%.8g PV dt=%.8g PV\n"
            "[SCW] independent mobile phase, IF97 density=%.12g kg/m3\n"
            "[INJECTION] reported q(H2O/CO2)=%.12g/%.12g mL/min; "
            "resolved reservoir total=%.12g m3/s\n",
            experiment.id, experiment.sourceExperiment,
            Sun2024Factorial::roleName(experiment.role),
            Sun2024Factorial::rateBasisName(study.rateBasis), study.nx,
            Sun2024Factorial::temperature, experiment.pressurePa,
            Sun2024Exp12::poreVolume, study.targetPoreVolumes,
            effectiveDtPv, waterDensity,
            experiment.waterRateMlPerMinute,
            experiment.co2RateMlPerMinute, totalRate);

    auto runtimeOptions =
        MPMC::cases::makeRuntimeOptions<Indices, Runtime, Config>(runOptions);
    runtimeOptions.scaling.enabled = true;
    runtimeOptions.scaling.rateWellResidualFloor = totalRate;
    runtimeOptions.scaling.validate();

    Runtime runtime(grid, fluid, runtimeOptions);
    runtime.setWells(makeWells(grid, definitions, injection));
    Vec solution = grid.createGlobalVector(Indices::numPrimaryVariables);
    initializeSolution(grid, solution, definitions, experiment.pressurePa);
    runtime.initializePhaseStateFromSolution(
        solution, MPMC::HydrocarbonPhaseState::LiquidOnly);
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
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscInitialize(&argc, &argv, nullptr, nullptr));
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
