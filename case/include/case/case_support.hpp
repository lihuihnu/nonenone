/**
 * @file case_support.hpp
 * @brief 算例公共运行框架：参数读取、时间推进、输出和诊断。
 */
#pragma once

#include <adaptive_timestep/adaptive_timestep.hpp>
#include <adaptive_timestep/natural/natural_petsc_backend.hpp>
#include <adaptive_timestep/well/control_cycle.hpp>
#include <case/case_validation.hpp>
#include <case/fluid_factory.hpp>
#include <case/failure_diagnostics.hpp>
#include <case/solver_diagnostics.hpp>
#include <common/console.hpp>
#include <common/petsc_io.hpp>
#include <common/units.hpp>
#include <natural/petsc/callbacks.hpp>
#include <output/simulation/model_inventory_output.hpp>
#include <output/simulation/component_mass_balance_output.hpp>
#include <output/simulation/reservoir_diagnostics.hpp>
#include <output/well/control_switch_output.hpp>
#include <output/well/detailed_well_output.hpp>
#include <output/well/producer_composition_output.hpp>

#include <petscksp.h>
#include <petscsnes.h>
#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace MPMC::cases
{

inline constexpr double secondsPerDay = 86400.0;
inline constexpr double gasConstant = MPMC::units::gasConstant; // J/(mol K)

/**
 * @brief 每个算例共用的少量运行参数。
 *
 * case_config.hpp 给默认值；命令行只负责临时覆盖，不复制物理参数。
 */
struct RunOptions final
{
    std::size_t numberOfSteps{1};
    std::size_t outputEvery{1};
    double dtDays{1.0};
    bool adaptive{true};
    double targetPVI{std::numeric_limits<double>::quiet_NaN()};
    std::string resultDirectory{"./results"};
    std::string meshDirectory;
};

/** @brief 探测算例是否提供局部 SNES 停滞配置，旧算例无需增加字段。 */
template <class Config, class = void>
struct HasSnesStagnationSettings : std::false_type
{};

template <class Config>
struct HasSnesStagnationSettings<
    Config,
    std::void_t<
        decltype(Config::Numerics::enableSnesStagnationGuard),
        decltype(Config::Numerics::snesStagnationMinimumIterations),
        decltype(Config::Numerics::snesStagnationWindow),
        decltype(Config::Numerics::snesStagnationRelativeImprovement)>> : std::true_type
{};


template <class Time, class = void>
struct HasTargetPVI : std::false_type
{};

template <class Time>
struct HasTargetPVI<Time, std::void_t<decltype(Time::targetPVI)>> : std::true_type
{};

template <class Time, class = void>
struct HasMaximumInternalDt : std::false_type
{};

template <class Time>
struct HasMaximumInternalDt<
    Time,
    std::void_t<decltype(Time::maximumDtDays)>> : std::true_type
{};

template <class Output, class = void>
struct HasProducerCompositionOutput : std::false_type
{};

template <class Output>
struct HasProducerCompositionOutput<
    Output,
    std::void_t<decltype(Output::writeProducerComposition)>> : std::true_type
{};

template <class Output, class = void>
struct HasProducerEffectivePoreVolume : std::false_type
{};

template <class Output>
struct HasProducerEffectivePoreVolume<
    Output,
    std::void_t<decltype(Output::producerEffectivePoreVolumeM3)>> : std::true_type
{};

template <class Fluid, class = void>
struct HasExplicitProducerLightHeavyGroups : std::false_type
{};

template <class Fluid>
struct HasExplicitProducerLightHeavyGroups<
    Fluid,
    std::void_t<
        decltype(Fluid::producerLighteningLightComponents),
        decltype(Fluid::producerLighteningHeavyComponents)>> : std::true_type
{};

template <class Fluid, class = void>
struct HasSingleLightHeavyComponents : std::false_type
{};

template <class Fluid>
struct HasSingleLightHeavyComponents<
    Fluid,
    std::void_t<
        decltype(Fluid::lightComponent),
        decltype(Fluid::heavyComponent)>> : std::true_type
{};


template <class Fluid, class = void>
struct HasWaterComponentIndex : std::false_type
{};

template <class Fluid>
struct HasWaterComponentIndex<
    Fluid,
    std::void_t<decltype(Fluid::waterComponent)>> : std::true_type
{};

/** @brief 构造非线性平台早停配置，并允许 PETSc 命令行临时覆盖。 */
template <class Config>
MPMC::NonlinearStagnationConfig makeNonlinearStagnationConfig()
{
    MPMC::NonlinearStagnationConfig cfg;
    if constexpr (HasSnesStagnationSettings<Config>::value)
    {
        cfg.enabled = Config::Numerics::enableSnesStagnationGuard;
        cfg.minimumIterations = Config::Numerics::snesStagnationMinimumIterations;
        cfg.stagnantIterations = Config::Numerics::snesStagnationWindow;
        cfg.requiredRelativeImprovement =
            Config::Numerics::snesStagnationRelativeImprovement;
    }

    PetscBool enabled = cfg.enabled ? PETSC_TRUE : PETSC_FALSE;
    PetscInt minimumIterations = static_cast<PetscInt>(cfg.minimumIterations);
    PetscInt stagnantIterations = static_cast<PetscInt>(cfg.stagnantIterations);
    PetscReal relativeImprovement =
        static_cast<PetscReal>(cfg.requiredRelativeImprovement);

    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetBool(
        nullptr, nullptr, "-snes_stagnation_guard", &enabled, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetInt(
        nullptr, nullptr, "-snes_stagnation_min_it", &minimumIterations, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetInt(
        nullptr, nullptr, "-snes_stagnation_window", &stagnantIterations, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsGetReal(
        nullptr, nullptr, "-snes_stagnation_rel_improvement", &relativeImprovement, nullptr));

    cfg.enabled = enabled == PETSC_TRUE;
    cfg.minimumIterations = static_cast<int>(minimumIterations);
    cfg.stagnantIterations = static_cast<int>(stagnantIterations);
    cfg.requiredRelativeImprovement = static_cast<double>(relativeImprovement);
    cfg.validate();
    return cfg;
}

template <class Config>
RunOptions readRunOptions()
{
    RunOptions options;
    options.numberOfSteps = Config::Time::numberOfSteps;
    options.outputEvery = Config::Output::every;
    options.dtDays = Config::Time::dtDays;
    options.adaptive = Config::Time::adaptive;
    if constexpr (HasTargetPVI<typename Config::Time>::value)
        options.targetPVI = Config::Time::targetPVI;
    options.resultDirectory = Config::Output::directory;

    // 所有算例都显式提供 Grid::meshDirectory：
    // StructuredGrid 使用空字符串，CpGrid 给出默认网格目录。
    // 运行配置使用 C++17 traits 表达，避免把编译器特性探测混入算例 API。
    options.meshDirectory = Config::Grid::meshDirectory;

    PetscInt steps = static_cast<PetscInt>(options.numberOfSteps);
    PetscInt outputEvery = static_cast<PetscInt>(options.outputEvery);
    PetscReal dtDays = options.dtDays;
    PetscBool adaptive = options.adaptive ? PETSC_TRUE : PETSC_FALSE;
    PetscReal targetPVI = options.targetPVI;
    char result[PETSC_MAX_PATH_LEN]{};
    char mesh[PETSC_MAX_PATH_LEN]{};
    PetscBool resultSet = PETSC_FALSE;
    PetscBool meshSet = PETSC_FALSE;

    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetInt(nullptr, nullptr, "-numSteps", &steps, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetInt(nullptr, nullptr, "-out_step", &outputEvery, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetReal(nullptr, nullptr, "-dt", &dtDays, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetBool(nullptr, nullptr, "-adaptive_dt", &adaptive, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetReal(nullptr, nullptr, "-target_pvi", &targetPVI, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetString(nullptr, nullptr, "-result_dir",
                                         result, sizeof(result), &resultSet));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetString(nullptr, nullptr, "-mesh_dir",
                                         mesh, sizeof(mesh), &meshSet));

    if (steps <= 0)
        throw std::invalid_argument("-numSteps must be positive.");
    if (outputEvery <= 0)
        throw std::invalid_argument("-out_step must be positive.");
    if (!(dtDays > 0.0) || !std::isfinite(dtDays))
        throw std::invalid_argument("-dt must be finite and positive (day).");

    options.numberOfSteps = static_cast<std::size_t>(steps);
    options.outputEvery = static_cast<std::size_t>(outputEvery);
    options.dtDays = dtDays;
    options.adaptive = adaptive == PETSC_TRUE;
    options.targetPVI = static_cast<double>(targetPVI);
    if (std::isfinite(options.targetPVI) && !(options.targetPVI > 0.0))
        throw std::invalid_argument("-target_pvi must be positive when supplied.");
    if (resultSet)
        options.resultDirectory = result;
    if (meshSet)
        options.meshDirectory = mesh;
    return options;
}

template <class Indices, class Config>
std::array<double, Indices::numComponents> makeStandardGasDensity()
{
    std::array<double, Indices::numComponents> result{};
    for (int c = 0; c < Indices::numComponents; ++c)
    {
        const std::size_t i = static_cast<std::size_t>(c);
        result[i] = Config::Fluid::molarMass[i] *
                    Config::Adsorption::standardPressure /
                    (gasConstant * Config::Adsorption::standardTemperature);
    }
    return result;
}

template <class Indices, class Runtime, class Config>
typename Runtime::Options makeRuntimeOptions(const RunOptions &run)
{
    typename Runtime::Options options;
    options.timeStep = run.dtDays * secondsPerDay;
    options.fugacityScalingFactor = Config::Numerics::fugacityScalingFactor;
    options.useVariableBounds = Config::Numerics::useVariableBounds;
    options.landConstant = Config::Land::constant;
    options.rockDensity = Config::Adsorption::rockDensity;

    if constexpr (Indices::hasAqueousCO2Dissolution)
        options.dissolvedCO2Component = Config::Dissolution::component;

    if constexpr (Indices::hasAdsorption)
        options.standardGasDensity = makeStandardGasDensity<Indices, Config>();

    return options;
}

template <class Config>
MPMC::AdaptiveTimeStepConfig makeTimeStepConfig(const RunOptions &run)
{
    MPMC::AdaptiveTimeStepConfig config;
    config.fixedOutputDt = run.dtDays * secondsPerDay;
    if constexpr (HasMaximumInternalDt<typename Config::Time>::value)
        config.maximumDt = Config::Time::maximumDtDays * secondsPerDay;
    config.minimumDt = Config::Time::minimumDtDays * secondsPerDay;
    config.cutFactor = Config::Time::cutFactor;
    config.growthFactor = Config::Time::growthFactor;
    config.difficultShrinkFactor = Config::Time::difficultShrinkFactor;
    config.easyNonlinearIterations = Config::Time::easyNonlinearIterations;
    config.difficultNonlinearIterations = Config::Time::difficultNonlinearIterations;
    config.maximumRetries = Config::Time::maximumRetries;
    config.maximumWellControlIterations = Config::Time::maximumWellControlIterations;
    config.adaptive = run.adaptive;

    // case_config.hpp 是默认值唯一来源；以下 PETSc 选项仅用于临时覆盖。
    PetscReal maximumDtDays = config.maximumDt / secondsPerDay;
    PetscReal minimumDtDays = config.minimumDt / secondsPerDay;
    PetscReal cutFactor = config.cutFactor;
    PetscReal growthFactor = config.growthFactor;
    PetscReal difficultShrink = config.difficultShrinkFactor;
    PetscInt easy = config.easyNonlinearIterations;
    PetscInt difficult = config.difficultNonlinearIterations;
    PetscInt retries = config.maximumRetries;
    PetscInt wellIterations = config.maximumWellControlIterations;

    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetReal(nullptr, nullptr, "-dt_max", &maximumDtDays, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetReal(nullptr, nullptr, "-dt_min", &minimumDtDays, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetReal(nullptr, nullptr, "-dt_cut", &cutFactor, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetReal(nullptr, nullptr, "-dt_growth", &growthFactor, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetReal(nullptr, nullptr, "-dt_difficult_shrink", &difficultShrink, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetInt(nullptr, nullptr, "-dt_easy_snes", &easy, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetInt(nullptr, nullptr, "-dt_difficult_snes", &difficult, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetInt(nullptr, nullptr, "-dt_max_retries", &retries, nullptr));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscOptionsGetInt(nullptr, nullptr, "-dt_max_well_control_iterations",
                                      &wellIterations, nullptr));

    config.maximumDt = static_cast<double>(maximumDtDays) * secondsPerDay;
    config.minimumDt = static_cast<double>(minimumDtDays) * secondsPerDay;
    config.cutFactor = static_cast<double>(cutFactor);
    config.growthFactor = static_cast<double>(growthFactor);
    config.difficultShrinkFactor = static_cast<double>(difficultShrink);
    config.easyNonlinearIterations = static_cast<int>(easy);
    config.difficultNonlinearIterations = static_cast<int>(difficult);
    config.maximumRetries = static_cast<int>(retries);
    config.maximumWellControlIterations = static_cast<int>(wellIterations);
    config.validate();
    return config;
}

/** @brief RAII 管理 SNES、残差 Vec 与 Jacobian Mat。 */
template <class Runtime>
class NaturalSolver final
{
public:
    NaturalSolver(Runtime &runtime, DM dm)
        : runtime_(runtime)
    {
        residual_ = runtime.createResidualVector();
        jacobian_ = runtime.createJacobian();
        PetscCallAbort(PETSC_COMM_WORLD, SNESCreate(PETSC_COMM_WORLD, &snes_));
        PetscCallAbort(PETSC_COMM_WORLD, SNESSetDM(snes_, dm));
        PetscCallAbort(PETSC_COMM_WORLD,
                       MPMC::installNaturalCallbacks(snes_, runtime, residual_, jacobian_));
        PetscCallAbort(PETSC_COMM_WORLD, SNESSetFromOptions(snes_));
        installMeshNormalizedConvergenceIfRequested_();
        // SNESSetFromOptions may replace the line-search implementation and
        // clear its user checks.  Bind Natural hooks to the final selected
        // line search before setup/solve.
        PetscCallAbort(PETSC_COMM_WORLD,
                       MPMC::installNaturalLineSearchCallbacks(snes_, runtime));
        PetscCallAbort(PETSC_COMM_WORLD, SNESSetUp(snes_));
    }

    ~NaturalSolver()
    {
        if (snes_) SNESDestroy(&snes_);
        if (jacobian_) MatDestroy(&jacobian_);
        if (residual_) VecDestroy(&residual_);
    }

    NaturalSolver(const NaturalSolver &) = delete;
    NaturalSolver &operator=(const NaturalSolver &) = delete;

    [[nodiscard]] SNES snes() const noexcept { return snes_; }

private:
    struct MeshNormalizedConvergenceContext final
    {
        Runtime *runtime{nullptr};
        PetscReal rmsAbsoluteTolerance{0.0};
        PetscReal infinityAbsoluteTolerance{0.0};
        PetscReal globalSignedMassAbsoluteTolerance{0.0};
        PetscInt globalEquationCount{0};
    };

    static PetscErrorCode meshNormalizedConvergenceTest_(
        SNES snes,
        PetscInt iteration,
        PetscReal xNorm,
        PetscReal stepNorm,
        PetscReal functionNorm,
        SNESConvergedReason *reason,
        void *context)
    {
        PetscFunctionBeginUser;
        PetscCheck(
            context != nullptr,
            PetscObjectComm(reinterpret_cast<PetscObject>(snes)),
            PETSC_ERR_ARG_NULL,
            "Mesh-normalized SNES convergence context is null.");

        auto &cfg =
            *static_cast<MeshNormalizedConvergenceContext *>(context);

        SNESConvergedReason defaultReason = SNES_CONVERGED_ITERATING;
        PetscCall(SNESConvergedDefault(
            snes,
            iteration,
            xNorm,
            stepNorm,
            functionNorm,
            &defaultReason,
            nullptr));

        *reason = defaultReason;
        if (defaultReason > 0)
        {
            const PetscReal rms =
                cfg.globalEquationCount > 0
                    ? functionNorm /
                        std::sqrt(static_cast<PetscReal>(cfg.globalEquationCount))
                    : std::numeric_limits<PetscReal>::infinity();
            if (rms > cfg.rmsAbsoluteTolerance)
            {
                *reason = SNES_CONVERGED_ITERATING;
                PetscFunctionReturn(PETSC_SUCCESS);
            }

            Vec residual = nullptr;
            PetscCall(SNESGetFunction(snes, &residual, nullptr, nullptr));
            PetscCheck(
                residual != nullptr,
                PetscObjectComm(reinterpret_cast<PetscObject>(snes)),
                PETSC_ERR_ARG_NULL,
                "SNES residual is unavailable for mesh-normalized convergence gates.");
            PetscCheck(
                cfg.runtime != nullptr,
                PetscObjectComm(reinterpret_cast<PetscObject>(snes)),
                PETSC_ERR_ARG_NULL,
                "Natural runtime is unavailable for global mass-residual convergence gate.");

            PetscReal infinityNorm = 0.0;
            PetscCall(VecNorm(residual, NORM_INFINITY, &infinityNorm));
            if (infinityNorm > cfg.infinityAbsoluteTolerance)
            {
                *reason = SNES_CONVERGED_ITERATING;
                PetscFunctionReturn(PETSC_SUCCESS);
            }

            const auto massResidual =
                cfg.runtime->evaluateGlobalSignedMassResidual(residual);
            if (massResidual.maximumAbsolute() >
                cfg.globalSignedMassAbsoluteTolerance)
            {
                *reason = SNES_CONVERGED_ITERATING;
                PetscFunctionReturn(PETSC_SUCCESS);
            }
        }

        PetscFunctionReturn(PETSC_SUCCESS);
    }

    void installMeshNormalizedConvergenceIfRequested_()
    {
        PetscReal rmsTolerance = 0.0;
        PetscReal infinityTolerance = 0.0;
        PetscReal massSumTolerance = 0.0;
        PetscBool rmsSet = PETSC_FALSE;
        PetscBool infinitySet = PETSC_FALSE;
        PetscBool massSumSet = PETSC_FALSE;
        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetReal(
                nullptr, nullptr, "-snes_rms_atol",
                &rmsTolerance, &rmsSet));
        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetReal(
                nullptr, nullptr, "-snes_linf_atol",
                &infinityTolerance, &infinitySet));
        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetReal(
                nullptr, nullptr, "-snes_global_mass_atol",
                &massSumTolerance, &massSumSet));

        const int requestedCount =
            (rmsSet == PETSC_TRUE ? 1 : 0) +
            (infinitySet == PETSC_TRUE ? 1 : 0) +
            (massSumSet == PETSC_TRUE ? 1 : 0);
        if (requestedCount != 0 && requestedCount != 3)
        {
            throw std::invalid_argument(
                "-snes_rms_atol, -snes_linf_atol and -snes_global_mass_atol "
                "must be supplied together.");
        }
        if (requestedCount == 0)
            return;
        if (!(rmsTolerance > 0.0) || !std::isfinite(rmsTolerance) ||
            !(infinityTolerance > 0.0) || !std::isfinite(infinityTolerance) ||
            !(massSumTolerance > 0.0) || !std::isfinite(massSumTolerance))
        {
            throw std::invalid_argument(
                "Mesh-normalized SNES RMS/Linf/global-mass tolerances "
                "must be finite and positive.");
        }

        PetscInt globalRows = 0;
        PetscCallAbort(
            PETSC_COMM_WORLD,
            VecGetSize(residual_, &globalRows));
        if (globalRows <= 0)
            throw std::logic_error(
                "Mesh-normalized SNES convergence requires a non-empty residual vector.");

        PetscReal oldAtol = 0.0;
        PetscReal rtol = 0.0;
        PetscReal stol = 0.0;
        PetscInt maxIterations = 0;
        PetscInt maxFunctions = 0;
        PetscCallAbort(
            PETSC_COMM_WORLD,
            SNESGetTolerances(
                snes_,
                &oldAtol,
                &rtol,
                &stol,
                &maxIterations,
                &maxFunctions));

        const PetscReal l2Tolerance =
            rmsTolerance * std::sqrt(static_cast<PetscReal>(globalRows));

        PetscCallAbort(
            PETSC_COMM_WORLD,
            SNESSetTolerances(
                snes_,
                l2Tolerance,
                rtol,
                stol,
                maxIterations,
                maxFunctions));

        meshNormalizedContext_.runtime = &runtime_;
        meshNormalizedContext_.rmsAbsoluteTolerance = rmsTolerance;
        meshNormalizedContext_.infinityAbsoluteTolerance =
            infinityTolerance;
        meshNormalizedContext_.globalSignedMassAbsoluteTolerance =
            massSumTolerance;
        meshNormalizedContext_.globalEquationCount = globalRows;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            SNESSetConvergenceTest(
                snes_,
                meshNormalizedConvergenceTest_,
                &meshNormalizedContext_,
                nullptr));

        PetscPrintf(
            PETSC_COMM_WORLD,
            "[SNES][MESH-NORM] N=%lld RMS_ATOL=%.12e L2_ATOL=%.12e "
            "LINF_ATOL=%.12e GLOBAL_MASS_ATOL=%.12e kg/s\n",
            static_cast<long long>(globalRows),
            static_cast<double>(rmsTolerance),
            static_cast<double>(l2Tolerance),
            static_cast<double>(infinityTolerance),
            static_cast<double>(massSumTolerance));
    }

    Runtime &runtime_;
    SNES snes_{nullptr};
    Vec residual_{nullptr};
    Mat jacobian_{nullptr};
    MeshNormalizedConvergenceContext meshNormalizedContext_{};
};

/**
 * @brief 将任意组分名称转换为紧凑且稳定的 CSV 列名。
 *
 * Component names in case files may contain `+`, `-`, `/` and spaces.  CSV itself
 * can quote such names, but normalized identifiers are much easier to use from
 * pandas/MATLAB/R.  Non-alphanumeric characters are therefore replaced by `_`.
 */
inline std::string csvColumnToken(std::string name)
{
    for (char &ch : name)
    {
        const unsigned char value = static_cast<unsigned char>(ch);
        if (!std::isalnum(value) && ch != '_')
            ch = '_';
    }
    return name;
}

/** @brief 返回单元 Natural 主变量块的描述性列名。 */
template <class Indices, class Config>
std::vector<std::string> solutionCsvColumns()
{
    std::vector<std::string> columns;
    columns.reserve(static_cast<std::size_t>(Indices::numPrimaryVariables));
    columns.emplace_back("pressure_Pa");

    auto appendCompositionColumns = [&](const char *prefix) {
        for (int component = 0;
             component < Indices::numIndependentCompositionsPerPhase;
             ++component)
        {
            columns.emplace_back(
                std::string(prefix) + std::to_string(component) + "_" +
                csvColumnToken(Config::Fluid::componentNames[
                    static_cast<std::size_t>(component)]));
        }
    };

    if constexpr (Indices::fullyCompositionalThreePhase)
    {
        appendCompositionColumns("oil_x_");
        appendCompositionColumns("gas_y_");
        appendCompositionColumns("water_x_");
        columns.emplace_back("oil_saturation");
        columns.emplace_back("gas_saturation");
        columns.emplace_back("water_saturation");
        if constexpr (Indices::hasWellUnknown)
            columns.emplace_back("well_pressure_Pa");
    }
    else
    {
        appendCompositionColumns("liquid_x_");

        if constexpr (Indices::hasWater)
            columns.emplace_back("water_saturation");
        if constexpr (Indices::hasWellUnknown)
            columns.emplace_back("well_pressure_Pa");

        columns.emplace_back("liquid_saturation");
        appendCompositionColumns("vapor_y_");
        columns.emplace_back("vapor_saturation");

        if constexpr (Indices::hasAqueousCO2Dissolution)
            columns.emplace_back("aqueous_CO2_mole_fraction");
    }

    if (columns.size() != static_cast<std::size_t>(Indices::numPrimaryVariables))
        throw std::logic_error("Solution CSV column layout does not match Indices.");
    return columns;
}

/** @brief 返回单元热力学相态块的描述性列名。 */
template <class Indices, class Config>
std::vector<std::string> phaseStateCsvColumns()
{
    std::vector<std::string> columns;
    columns.reserve(static_cast<std::size_t>(Indices::numPhaseStateVariables));

    if constexpr (Indices::fullyCompositionalThreePhase)
    {
        columns.emplace_back("phase_presence_mask");

        for (int component = 0; component < Indices::numComponents; ++component)
        {
            columns.emplace_back(
                "K_gas_over_oil_" + std::to_string(component) + "_" +
                csvColumnToken(Config::Fluid::componentNames[
                    static_cast<std::size_t>(component)]));
        }
        for (int component = 0; component < Indices::numComponents; ++component)
        {
            columns.emplace_back(
                "K_water_over_oil_" + std::to_string(component) + "_" +
                csvColumnToken(Config::Fluid::componentNames[
                    static_cast<std::size_t>(component)]));
        }
        for (int component = 0; component < Indices::numComponents; ++component)
        {
            columns.emplace_back(
                "z_" + std::to_string(component) + "_" +
                csvColumnToken(Config::Fluid::componentNames[
                    static_cast<std::size_t>(component)]));
        }

        columns.emplace_back("beta_oil");
        columns.emplace_back("beta_gas");
        columns.emplace_back("beta_water");
        columns.emplace_back("Z_oil");
        columns.emplace_back("Z_gas");
        columns.emplace_back("Z_water");
        columns.emplace_back("phase_hysteresis_suppression_mask");
    }
    else
    {
        columns.emplace_back("phase_flag");

        for (int component = 0; component < Indices::numComponents; ++component)
        {
            columns.emplace_back(
                "K_" + std::to_string(component) + "_" +
                csvColumnToken(Config::Fluid::componentNames[
                    static_cast<std::size_t>(component)]));
        }

        for (int component = 0; component < Indices::numComponents; ++component)
        {
            columns.emplace_back(
                "z_" + std::to_string(component) + "_" +
                csvColumnToken(Config::Fluid::componentNames[
                    static_cast<std::size_t>(component)]));
        }

        columns.emplace_back("liquid_mole_fraction");
        columns.emplace_back("liquid_compressibility_factor");
        columns.emplace_back("vapor_compressibility_factor");
    }

    if (columns.size() != static_cast<std::size_t>(Indices::numPhaseStateVariables))
        throw std::logic_error("Phase-state CSV column layout does not match Indices.");
    return columns;
}

/**
 * @brief 按网格原始输入顺序将单元 runtime Vec 保存为 CSV。
 *
 * Natural/CpGrid vectors are intentionally stored internally by partition/current id
 * so every MPI rank owns a contiguous PETSc range.  External result files must not
 * expose that partition-dependent ordering.  The runtime therefore creates a
 * temporary input-ordered copy immediately before CSV output.
 *
 * CSV layout is one row per input cell.  Column 0 is `input_index`; the remaining
 * columns are the DOFs of that cell in the exact Natural layout supplied by
 * @p columnNames.
 */
template <class Runtime>
void saveInputOrderedCsv(
    Runtime &runtime,
    Vec currentOrderedVector,
    PetscInt dofPerCell,
    const std::string &filename,
    const std::vector<std::string> &columnNames)
{
    Vec inputOrdered = runtime.createInputOrderedCopy(
        currentOrderedVector,
        dofPerCell);

    PetscCallAbort(
        PETSC_COMM_WORLD,
        MPMC::petsc::saveCellVectorCsv(
            inputOrdered,
            dofPerCell,
            filename,
            columnNames));

    PetscCallAbort(
        PETSC_COMM_WORLD,
        VecDestroy(&inputOrdered));
}

/** Print one uniform model/time/solver report immediately before time stepping. */
template <class Indices, class Runtime, class Config>
void printSimulationConfiguration(
    const Runtime &runtime,
    SNES snes,
    const RunOptions &run,
    const MPMC::AdaptiveTimeStepConfig &timeConfig,
    const MPMC::NonlinearStagnationConfig &stagnationConfig,
    PetscMPIInt rank)
{
    if (rank != 0)
        return;

    const auto &fluid = runtime.fluidSystem();
    std::ostringstream components;
    for (int c = 0; c < Indices::numComponents; ++c)
    {
        if (c != 0)
            components << ", ";
        components << fluid.componentNames[static_cast<std::size_t>(c)];
    }

    const char *snesType = "unknown";
    PetscReal snesAtol = 0.0;
    PetscReal snesRtol = 0.0;
    PetscReal snesStol = 0.0;
    PetscInt snesMaxIts = 0;
    PetscInt snesMaxFuncs = 0;
    PetscCallAbort(PETSC_COMM_WORLD, SNESGetType(snes, &snesType));
    PetscCallAbort(PETSC_COMM_WORLD, SNESGetTolerances(
        snes, &snesAtol, &snesRtol, &snesStol, &snesMaxIts, &snesMaxFuncs));

    KSP ksp = nullptr;
    PetscCallAbort(PETSC_COMM_WORLD, SNESGetKSP(snes, &ksp));
    const char *kspType = "unknown";
    PetscReal kspRtol = 0.0;
    PetscReal kspAtol = 0.0;
    PetscReal kspDtol = 0.0;
    PetscInt kspMaxIts = 0;
    PC pc = nullptr;
    const char *pcType = "unknown";
    if (ksp != nullptr)
    {
        PetscCallAbort(PETSC_COMM_WORLD, KSPGetType(ksp, &kspType));
        PetscCallAbort(PETSC_COMM_WORLD, KSPGetTolerances(
            ksp, &kspRtol, &kspAtol, &kspDtol, &kspMaxIts));
        PetscCallAbort(PETSC_COMM_WORLD, KSPGetPC(ksp, &pc));
        if (pc != nullptr)
            PetscCallAbort(PETSC_COMM_WORLD, PCGetType(pc, &pcType));
    }

    MPMC::ConsoleSection section("SIMULATION CONFIGURATION");
    section.row("Case", Config::name)
        .row("Cells", std::to_string(runtime.cellCount()))
        .row("Components / phases", std::to_string(Indices::numComponents) + " / " +
            std::to_string(Indices::numPhases))
        .row("Primary unknowns / cell", std::to_string(Indices::numPrimaryVariables))
        .row("Phase formulation", Indices::fullyCompositionalThreePhase
            ? "fully compositional O/G/W"
            : "legacy oil-gas + independent water")
        .row("Thermodynamic backend",
             MPMC::thermodynamicModelName(fluid.eos.thermodynamicModel()))
        .row("Temperature", MPMC::consoleNumber(fluid.temperature, 3), "K")
        .row("Components", components.str())
        .row("Wells", std::to_string(runtime.wells().size()))
        .row("Aqueous CO2 dissolution", MPMC::consoleOnOff(Indices::hasAqueousCO2Dissolution))
        .row("Adsorption", MPMC::consoleOnOff(Indices::hasAdsorption))
        .row("Land trapping", MPMC::consoleOnOff(Indices::hasLandTrapping))
        .row("Fugacity residual scale", MPMC::consoleNumber(Config::Numerics::fugacityScalingFactor))
        .row("Variable bounds", MPMC::consoleOnOff(Config::Numerics::useVariableBounds))
        .separator()
        .row("Output intervals", std::to_string(run.numberOfSteps))
        .row("Output interval", MPMC::consoleNumber(timeConfig.fixedOutputDt / secondsPerDay), "day")
        .row("Adaptive time stepping", MPMC::consoleOnOff(timeConfig.adaptive))
        .row("Maximum internal dt", MPMC::consoleNumber(timeConfig.effectiveMaximumDt() / secondsPerDay), "day")
        .row("Minimum dt", MPMC::consoleNumber(timeConfig.minimumDt / secondsPerDay), "day")
        .row("dt cut / growth", MPMC::consoleNumber(timeConfig.cutFactor, 4) + " / " +
            MPMC::consoleNumber(timeConfig.growthFactor, 4))
        .row("Easy / difficult SNES",
             std::to_string(timeConfig.easyNonlinearIterations) + " / " +
             std::to_string(timeConfig.difficultNonlinearIterations))
        .row("Maximum dt retries", std::to_string(timeConfig.maximumRetries))
        .row("Maximum well re-solves", std::to_string(timeConfig.maximumWellControlIterations))
        .separator()
        .row("SNES", snesType)
        .row("SNES atol / rtol", MPMC::consoleScientific(snesAtol, 3) + " / " +
            MPMC::consoleScientific(snesRtol, 3))
        .row("SNES max iterations", std::to_string(snesMaxIts))
        .row("SNES stagnation guard", MPMC::consoleOnOff(stagnationConfig.enabled))
        .row("Stagnation min-it / window",
             std::to_string(stagnationConfig.minimumIterations) + " / " +
             std::to_string(stagnationConfig.stagnantIterations))
        .row("Stagnation relative improvement",
             MPMC::consoleScientific(stagnationConfig.requiredRelativeImprovement, 3))
        .row("KSP / PC", std::string(kspType) + " / " + pcType)
        .row("KSP atol / rtol", MPMC::consoleScientific(kspAtol, 3) + " / " +
            MPMC::consoleScientific(kspRtol, 3))
        .row("KSP max iterations", std::to_string(kspMaxIts))
        .row("Results directory", run.resultDirectory);

    const std::string text = section.str();
    PetscPrintf(PETSC_COMM_SELF, "%s", text.c_str());
}

/** @brief 所有算例共用的收敛状态输出入口。 */
template <class Indices, class Runtime, class Config>
class CaseOutput final
{
public:
    CaseOutput(const RunOptions &run, PetscMPIInt rank)
        : directory_(run.resultDirectory),
          snapshotEvery_(run.outputEvery),
          rank_(rank),
          well_(wellOptions_(run), rank),
          producer_(producerOptions_(run), rank),
          inventory_(inventoryOptions_(run), rank),
          massBalance_(massBalanceOptions_(run), rank),
          reservoir_(reservoirOptions_(run), rank)
    {
        // 并行：仅 rank 0 创建目录/文件；进入 Vec collective I/O 前先同步。
        MPI_Barrier(PETSC_COMM_WORLD);
    }

    void write(std::size_t step, double time, Runtime &runtime, Vec solution)
    {
        well_.write(step, time, runtime, solution);
        producer_.write(step, time, runtime, solution);

        // 性能：组分守恒在每个固定输出状态已经计算全局库存，这里直接复用该
        // reduction，避免对同一接受解重复构造 EOS/单元物性并再次 MPI 归约。
        if constexpr (Config::Output::enableComponentMassBalance)
        {
            const auto reducedInventory =
                runtime.evaluateGlobalInventory(solution);
            inventory_.writeInventory(step, time, reducedInventory);
            massBalance_.writeInventory(step, time, reducedInventory);
        }
        else
        {
            inventory_.write(step, time, runtime, solution);
            massBalance_.write(step, time, runtime, solution);
        }

        reservoir_.write(step, time, runtime, solution);
        writeThermodynamicProfile_(step, time, runtime);

        if (!Config::Output::writeSolutionSnapshots)
            return;
        if (step != 0 && step % snapshotEvery_ != 0)
            return;

        const std::string suffix = "_step_" + std::to_string(step) + ".csv";
        saveInputOrderedCsv(
            runtime,
            solution,
            static_cast<PetscInt>(Indices::numPrimaryVariables),
            directory_ + "/solution" + suffix,
            solutionCsvColumns<Indices, Config>());

        if constexpr (Config::Output::writePhaseStateSnapshots)
        {
            saveInputOrderedCsv(
                runtime,
                runtime.phaseStateVector(),
                static_cast<PetscInt>(Indices::numPhaseStateVariables),
                directory_ + "/phase_state" + suffix,
                phaseStateCsvColumns<Indices, Config>());
        }
    }

    void acceptedStep(Runtime &runtime, Vec solution, double acceptedTime)
    {
        producer_.acceptedStep(runtime, solution, acceptedTime);
        massBalance_.acceptedStep(runtime, solution, acceptedTime);
    }

    [[nodiscard]] double producerPVI() const noexcept
    {
        return producer_.pvi();
    }

private:
    void writeThermodynamicProfile_(std::size_t step, double time, Runtime &runtime)
    {
        const auto &eos = runtime.fluidSystem().eos;
        if (!eos.thermodynamicProfilerEnabled() || !eos.usesCubicPlusAssociation())
            return;

        const auto local = eos.thermodynamicProfile();
        std::array<unsigned long long, 10> localCounts{
            local.cpaPhaseResultCalls,
            local.cpaMixingCalls,
            local.cpaTemperatureCacheHits,
            local.cpaAssociationCalls,
            local.cpaAssociationAnalyticCalls,
            local.cpaAssociationIterativeCalls,
            local.cpaAssociationIterations,
            local.cpaDensityRootCalls,
            local.cpaDensityResidualEvaluations,
            local.cpaDensityBisectionIterations};
        std::array<unsigned long long, 10> sumCounts{};
        PetscCallMPIAbort(
            PETSC_COMM_WORLD,
            MPI_Allreduce(
                localCounts.data(), sumCounts.data(),
                static_cast<int>(localCounts.size()),
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, PETSC_COMM_WORLD));

        std::array<double, 4> localSeconds{
            local.cpaPhaseResultSeconds,
            local.cpaMixingSeconds,
            local.cpaAssociationSeconds,
            local.cpaDensityRootSeconds};
        std::array<double, 4> maxSeconds{};
        std::array<double, 4> sumSeconds{};
        PetscCallMPIAbort(
            PETSC_COMM_WORLD,
            MPI_Allreduce(
                localSeconds.data(), maxSeconds.data(),
                static_cast<int>(localSeconds.size()),
                MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD));
        PetscCallMPIAbort(
            PETSC_COMM_WORLD,
            MPI_Allreduce(
                localSeconds.data(), sumSeconds.data(),
                static_cast<int>(localSeconds.size()),
                MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD));

        if (rank_ != 0)
            return;

        const double cacheHitRatio = sumCounts[1] > 0
            ? static_cast<double>(sumCounts[2]) / static_cast<double>(sumCounts[1])
            : 0.0;
        MPMC::ConsoleSection section("THERMODYNAMIC PROFILER");
        section.row("output step / time", std::to_string(step) + " / " +
                    MPMC::consoleNumber(time / secondsPerDay, 9) + " day")
            .row("CPA phase evaluations (sum ranks)", std::to_string(sumCounts[0]))
            .row("CPA mixing calls / cache hits",
                 std::to_string(sumCounts[1]) + " / " + std::to_string(sumCounts[2]) +
                 " (" + MPMC::consoleNumber(100.0 * cacheHitRatio, 2) + "%)")
            .row("Association calls analytic/iterative",
                 std::to_string(sumCounts[3]) + "  " +
                 std::to_string(sumCounts[4]) + " / " + std::to_string(sumCounts[5]))
            .row("Association site-solver iterations", std::to_string(sumCounts[6]))
            .row("Density roots / residual evals / bisections",
                 std::to_string(sumCounts[7]) + " / " +
                 std::to_string(sumCounts[8]) + " / " + std::to_string(sumCounts[9]))
            .separator()
            .row("max-rank CPA phase wall", MPMC::consoleNumber(maxSeconds[0], 6), "s")
            .row("max-rank mixing wall", MPMC::consoleNumber(maxSeconds[1], 6), "s")
            .row("max-rank association wall", MPMC::consoleNumber(maxSeconds[2], 6), "s")
            .row("max-rank density-root wall", MPMC::consoleNumber(maxSeconds[3], 6), "s")
            .row("sum-rank CPA phase wall", MPMC::consoleNumber(sumSeconds[0], 6), "s");
        const std::string text = section.str();
        PetscPrintf(PETSC_COMM_SELF, "%s", text.c_str());

        const auto path = std::filesystem::path(directory_) / "thermodynamic_profile.csv";
        const bool writeHeader = !std::filesystem::exists(path);
        std::ofstream out(path, std::ios::app);
        if (!out)
            throw std::runtime_error("Cannot open thermodynamic_profile.csv for writing.");
        if (writeHeader)
        {
            out << "output_step,time_day,phase_calls,mixing_calls,cache_hits,"
                   "association_calls,association_analytic_calls,association_iterative_calls,"
                   "association_iterations,density_root_calls,density_residual_evals,"
                   "density_bisection_iterations,phase_s_max,mixing_s_max,association_s_max,"
                   "density_root_s_max,phase_s_sum,mixing_s_sum,association_s_sum,density_root_s_sum\n";
        }
        out << step << ',' << std::setprecision(17) << time / secondsPerDay;
        for (auto value : sumCounts) out << ',' << value;
        for (auto value : maxSeconds) out << ',' << value;
        for (auto value : sumSeconds) out << ',' << value;
        out << '\n';
    }

    static MPMC::DetailedWellOutputOptions wellOptions_(const RunOptions &run)
    {
        MPMC::DetailedWellOutputOptions options;
        options.resultDirectory = run.resultDirectory;
        options.printState = Config::Output::printWells;
        options.printPhaseDetails = Config::Output::printWellPhaseDetails;
        options.writeHistory = Config::Output::writeWellHistory;
        options.printEveryOutputSteps = run.outputEvery;
        options.componentNames.reserve(
            static_cast<std::size_t>(Indices::numComponents));
        for (int component = 0;
             component < Indices::numComponents; ++component)
        {
            options.componentNames.emplace_back(
                Config::Fluid::componentNames[
                    static_cast<std::size_t>(component)]);
        }
        return options;
    }

    static MPMC::ProducerCompositionOutputOptions producerOptions_(
        const RunOptions &run)
    {
        MPMC::ProducerCompositionOutputOptions options;
        options.resultDirectory = run.resultDirectory;
        if constexpr (HasProducerCompositionOutput<
                          typename Config::Output>::value)
        {
            options.enabled =
                Config::Output::writeProducerComposition;
        }

        options.componentNames.reserve(
            static_cast<std::size_t>(Indices::numComponents));
        for (int component = 0;
             component < Indices::numComponents;
             ++component)
        {
            options.componentNames.emplace_back(
                Config::Fluid::componentNames[
                    static_cast<std::size_t>(component)]);
        }

        if constexpr (HasExplicitProducerLightHeavyGroups<
                          typename Config::Fluid>::value)
        {
            options.lightComponents.assign(
                Config::Fluid::producerLighteningLightComponents.begin(),
                Config::Fluid::producerLighteningLightComponents.end());
            options.heavyComponents.assign(
                Config::Fluid::producerLighteningHeavyComponents.begin(),
                Config::Fluid::producerLighteningHeavyComponents.end());
        }
        else if constexpr (HasSingleLightHeavyComponents<
                               typename Config::Fluid>::value)
        {
            options.lightComponents = {
                Config::Fluid::lightComponent};
            options.heavyComponents = {
                Config::Fluid::heavyComponent};
        }

        if constexpr (HasWaterComponentIndex<
                          typename Config::Fluid>::value)
        {
            options.waterComponent =
                Config::Fluid::waterComponent;
        }

        if constexpr (HasProducerEffectivePoreVolume<
                          typename Config::Output>::value)
        {
            options.effectivePoreVolumeM3 =
                Config::Output::producerEffectivePoreVolumeM3;
        }
        return options;
    }

    static MPMC::ModelInventoryOutputOptions inventoryOptions_(const RunOptions &run)
    {
        MPMC::ModelInventoryOutputOptions options;
        options.resultDirectory = run.resultDirectory;
        options.printState = Config::Output::printInventory;
        options.writeHistory = Config::Output::writeInventoryHistory;
        options.writeMassTotals = Config::Output::writeMassTotals;
        options.printEveryOutputSteps = run.outputEvery;
        if constexpr (Indices::hasAqueousCO2Dissolution)
            options.dissolvedCO2Component = Config::Dissolution::component;
        return options;
    }

    static MPMC::ComponentMassBalanceOutputOptions massBalanceOptions_(const RunOptions &run)
    {
        MPMC::ComponentMassBalanceOutputOptions options;
        options.resultDirectory = run.resultDirectory;
        options.enabled = Config::Output::enableComponentMassBalance;
        options.printState = Config::Output::printComponentMassBalance;
        options.writeHistory = Config::Output::writeComponentMassBalance;
        options.printEveryOutputSteps = run.outputEvery;

        PetscBool auditInternalFaces = PETSC_FALSE;
        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetBool(
                nullptr, nullptr,
                "-audit_internal_face_conservation",
                &auditInternalFaces, nullptr));
        options.auditInternalFaceConservation =
            auditInternalFaces == PETSC_TRUE;

        options.componentNames.reserve(static_cast<std::size_t>(Indices::numComponents));
        for (int c = 0; c < Indices::numComponents; ++c)
            options.componentNames.emplace_back(
                Config::Fluid::componentNames[static_cast<std::size_t>(c)]);
        return options;
    }

    static MPMC::ReservoirDiagnosticsOutputOptions reservoirOptions_(const RunOptions &run)
    {
        MPMC::ReservoirDiagnosticsOutputOptions options;
        options.resultDirectory = run.resultDirectory;
        options.printState = Config::Output::printReservoirDiagnostics;
        options.writeHistory = Config::Output::writeReservoirDiagnostics;
        options.printEveryOutputSteps = run.outputEvery;
        return options;
    }

    std::string directory_;
    std::size_t snapshotEvery_{1};
    PetscMPIInt rank_{0};
    MPMC::DetailedWellOutput<Indices, Runtime> well_;
    MPMC::ProducerCompositionOutput<Indices, Runtime> producer_;
    MPMC::ModelInventoryOutput<Indices, Runtime> inventory_;
    MPMC::ComponentMassBalanceOutput<Indices, Runtime> massBalance_;
    MPMC::ReservoirDiagnosticsOutput<Indices, Runtime> reservoir_;
};

template <class Config>
MPMC::WellControlSwitchOutputOptions makeWellSwitchOutputOptions(const RunOptions &run)
{
    MPMC::WellControlSwitchOutputOptions options;
    options.resultDirectory = run.resultDirectory;
    options.print = Config::Output::printWellControlSwitches;
    options.writeHistory = Config::Output::writeWellControlSwitchHistory;
    return options;
}

/**
 * @brief 所有算例通过统一的 AdaptiveTimeStepper 与诊断流程推进时间。
 */
template <class Indices, class Runtime, class Config>
void runTimeLoop(Runtime &runtime, SNES snes, Vec solution,
                 const RunOptions &run, PetscMPIInt rank)
{
    const auto timeConfig = makeTimeStepConfig<Config>(run);
    const auto stagnationConfig = makeNonlinearStagnationConfig<Config>();
    printSimulationConfiguration<Indices, Runtime, Config>(
        runtime, snes, run, timeConfig, stagnationConfig, rank);

    CaseOutput<Indices, Runtime, Config> output(run, rank);
    MPMC::WellControlSwitchOutput<Indices> switchOutput(
        makeWellSwitchOutputOptions<Config>(run), rank);

    double currentTime = 0.0;
    runtime.setCurrentTime(currentTime);
    output.write(0, currentTime, runtime, solution);

    using Well = typename Runtime::Well;
    auto stateProvider = [&runtime, solution](std::vector<MPMC::WellState<Indices>> &states) {
        states = runtime.evaluateWellStates(solution);
    };
    auto sync = [](const std::vector<Well> &) {};
    auto switchCallback = [&runtime, &switchOutput](
        const Well &well,
        const MPMC::WellState<Indices> &state,
        const MPMC::WellControlUpdate &update) {
        switchOutput.write(runtime.currentTime(), well, state, update);
    };

    using ControlCycle = MPMC::WellControlCycle<
        Indices,
        Well,
        decltype(stateProvider),
        decltype(sync),
        decltype(switchCallback)>;
    ControlCycle controls(
        runtime.mutableWells(), stateProvider, sync, switchCallback);

    auto acceptedStepHook = [&output](Runtime &acceptedRuntime, Vec acceptedSolution,
                                      double acceptedTime) {
        output.acceptedStep(acceptedRuntime, acceptedSolution, acceptedTime);
    };

    using FailedSolveHook = FailedSolveDiagnostics<Indices, Runtime, Config>;
    using Backend = MPMC::NaturalAdaptiveBackend<
        Runtime, ControlCycle, decltype(acceptedStepHook), FailedSolveHook>;
    using Observer = SolverDiagnosticsObserver;

    Backend backend(
        runtime,
        snes,
        solution,
        currentTime,
        std::move(controls),
        acceptedStepHook,
        FailedSolveHook(
            rank,
            run.resultDirectory,
            Config::Time::printAdaptiveSteps,
            Config::Output::writeSolverHistory),
        Config::Output::printNewtonIterations,
        stagnationConfig);
    MPMC::AdaptiveTimeStepper<Backend, Observer> stepper(
        backend,
        timeConfig,
        Observer(
            rank,
            SolverDiagnosticsOptions{
                run.resultDirectory,
                Config::Time::printAdaptiveSteps,
                Config::Output::writeSolverHistory,
                secondsPerDay,
                "day"}));

    PetscLogDouble localStart = 0.0;
    PetscLogDouble localEnd = 0.0;
    PetscCallAbort(PETSC_COMM_WORLD, PetscTime(&localStart));

    std::size_t completedFixedSteps = 0;
    if (std::isfinite(run.targetPVI))
    {
        completedFixedSteps = stepper.runUntil(
            run.numberOfSteps,
            [&output, targetPVI = run.targetPVI](std::size_t, double) {
                const double pvi = output.producerPVI();
                return std::isfinite(pvi) && pvi >= targetPVI;
            },
            [&output, &runtime, solution](std::size_t step, double time) {
                output.write(step, time, runtime, solution);
            });
    }
    else
    {
        stepper.run(
            run.numberOfSteps,
            [&output, &runtime, solution](std::size_t step, double time) {
                output.write(step, time, runtime, solution);
            });
        completedFixedSteps = run.numberOfSteps;
    }

    saveInputOrderedCsv(
        runtime,
        solution,
        static_cast<PetscInt>(Indices::numPrimaryVariables),
        run.resultDirectory + "/solution_final.csv",
        solutionCsvColumns<Indices, Config>());

    if constexpr (Config::Output::writePhaseStateSnapshots)
    {
        saveInputOrderedCsv(
            runtime,
            runtime.phaseStateVector(),
            static_cast<PetscInt>(Indices::numPhaseStateVariables),
            run.resultDirectory + "/phase_state_final.csv",
            phaseStateCsvColumns<Indices, Config>());
    }

    PetscCallAbort(PETSC_COMM_WORLD, PetscTime(&localEnd));
    const double localWall = static_cast<double>(localEnd - localStart);
    double globalWall = localWall;
    PetscCallMPIAbort(
        PETSC_COMM_WORLD,
        MPI_Allreduce(&localWall, &globalWall, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD));

    writeSimulationSummary(
        SimulationSummaryInfo{
            Config::name,
            run.resultDirectory,
            completedFixedSteps,
            runtime.currentTime() / secondsPerDay,
            switchOutput.switchCount(),
            Config::Output::printFinalSummary,
            Config::Output::writeFinalSummary},
        stepper.statistics(),
        globalWall,
        rank);
}

inline void printRunSummary(const char *name, const RunOptions &run, PetscMPIInt rank)
{
    if (rank != 0)
        return;

    MPMC::ConsoleSection section("RUN REQUEST");
    section.row("Case", name)
        .row("Output intervals", std::to_string(run.numberOfSteps))
        .row("Nominal output dt", MPMC::consoleNumber(run.dtDays), "day")
        .row("Adaptive dt", MPMC::consoleOnOff(run.adaptive))
        .row("Console OUT_STEP", std::to_string(run.outputEvery))
        .row("Results", run.resultDirectory);
    if (!run.meshDirectory.empty())
        section.row("Grid input", run.meshDirectory);

    const std::string text = section.str();
    PetscPrintf(PETSC_COMM_SELF, "%s", text.c_str());
}

} // namespace MPMC::cases
