/**
 * @file output_feature_test.cpp
 * @brief 单元测试：验证 `output_feature` 的核心语义、边界条件和回归行为。
 */
#include <output/output.hpp>

#include <indices/indices.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace
{

template <class Indices>
void requireExists(
    const std::filesystem::path &directory,
    const char *name,
    bool expected)
{
    const bool exists =
        std::filesystem::exists(
            directory / name);

    if (exists != expected)
    {
        throw std::runtime_error(
            std::string("Unexpected output feature file state: ") +
            name);
    }
}

template <class Indices>
void runConfig(
    const std::filesystem::path &directory)
{
    std::filesystem::remove_all(
        directory);

    {
        MPMC::OutputWriter<Indices>
            writer(directory);
    }

    requireExists<Indices>(
        directory,
        "solution_samples.csv",
        true);

    requireExists<Indices>(
        directory,
        "phase_state_samples.csv",
        true);

    requireExists<Indices>(
        directory,
        "well_phase_rates.csv",
        true);

    requireExists<Indices>(
        directory,
        "adsorbed_mass_samples.csv",
        Indices::hasAdsorption);

    requireExists<Indices>(
        directory,
        "adsorbed_mass_totals.csv",
        Indices::hasAdsorption);

    requireExists<Indices>(
        directory,
        "trapped_mass_samples.csv",
        Indices::hasLandTrapping);

    requireExists<Indices>(
        directory,
        "trapped_mass_totals.csv",
        Indices::hasLandTrapping);

    requireExists<Indices>(
        directory,
        "dissolved_mass_samples.csv",
        Indices::hasAqueousCO2Dissolution);

    requireExists<Indices>(
        directory,
        "dissolved_mass_totals.csv",
        Indices::hasAqueousCO2Dissolution);

    requireExists<Indices>(
        directory,
        "oil_phase_component_mass_samples.csv",
        true);

    requireExists<Indices>(
        directory,
        "oil_phase_component_mass_totals.csv",
        true);

    std::filesystem::remove_all(
        directory);
}

using Base =
    MPMC::ScalarIndices<
        MPMC::CompositionalModelConfig<
            6, true, true,
            false, false, false>>;

using Dissolution =
    MPMC::ScalarIndices<
        MPMC::CompositionalModelConfig<
            6, true, true,
            true, false, false>>;

using Land =
    MPMC::ScalarIndices<
        MPMC::CompositionalModelConfig<
            6, true, true,
            false, false, true>>;

using Adsorption =
    MPMC::ScalarIndices<
        MPMC::CompositionalModelConfig<
            6, true, true,
            false, true, false>>;

using Full =
    MPMC::ScalarIndices<
        MPMC::CompositionalModelConfig<
            6, true, true,
            true, true, true>>;

} // namespace

int main()
{
    const auto base =
        std::filesystem::temp_directory_path();

    runConfig<Base>(
        base / "mpmc_output_base");

    runConfig<Dissolution>(
        base / "mpmc_output_diss");

    runConfig<Land>(
        base / "mpmc_output_land");

    runConfig<Adsorption>(
        base / "mpmc_output_ads");

    runConfig<Full>(
        base / "mpmc_output_full");

    std::cout
        << "Output feature matrix: ALL PASS\n";

    return 0;
}
