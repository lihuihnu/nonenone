/**
 * @file output_test.cpp
 * @brief 单元测试：验证 `output` 的核心语义、边界条件和回归行为。
 */
#include <output/output.hpp>

#include <indices/indices.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using Config =
    MPMC::CompositionalModelConfig<
        6,
        true,
        true,
        true,
        true,
        true>;

using Indices =
    MPMC::ScalarIndices<Config>;

std::string readAll(
    const std::filesystem::path &path)
{
    std::ifstream file(path);
    if (!file)
        throw std::runtime_error(
            "Cannot read test output file.");

    std::ostringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

void requireEqual(
    const std::string &actual,
    const std::string &expected,
    const char *label)
{
    if (actual != expected)
    {
        std::cerr
            << "[FAIL] "
            << label
            << "\nEXPECTED:\n"
            << expected
            << "\nACTUAL:\n"
            << actual
            << '\n';

        throw std::runtime_error(
            "Output text mismatch.");
    }
}

} // namespace

int main()
{
    const MPMC::OutputCellSelection defaultCell{};
    if (defaultCell.space != MPMC::OutputCellIdSpace::InputIndex ||
        MPMC::OutputCellIdSpace::StorageIndex != MPMC::OutputCellIdSpace::InputIndex)
    {
        throw std::runtime_error(
            "Output cell numbering must default to canonical input-index space.");
    }

    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() /
        "mpmc_output_v2_core_test";

    std::filesystem::remove_all(directory);

    {
        MPMC::OutputWriter<Indices>
            writer(directory);

        writer.writeSolutionSamples(
            {1.25, -2.5});

        writer.writePhaseStateSamples(
            {0.125});

        writer.writeWellPhaseRates(
            std::array<double, Indices::numPhases>{
                -2.0,
                -200.0,
                -6.0});

        MPMC::ComponentTotals<6>
            totals;

        totals.component =
            {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

        writer.writeMassTotals(
            MPMC::MassSeries::Adsorbed,
            12.5,
            totals);

        writer.writeMassTotals(
            MPMC::MassSeries::Trapped,
            12.5,
            totals);

        writer.writeMassTotals(
            MPMC::MassSeries::Dissolved,
            12.5,
            totals);

        writer.writeMassTotals(
            MPMC::MassSeries::OilPhaseComponents,
            12.5,
            totals);

        writer.writeMassSamples(
            MPMC::MassSeries::Trapped,
            {1.0e-6, 2.0e-6});
    }

    requireEqual(
        readAll(
            directory / "solution_samples.csv"),
        "value_0,value_1\n"
        "1.25,-2.5\n",
        "solution_samples.csv");

    requireEqual(
        readAll(
            directory / "phase_state_samples.csv"),
        "value_0\n"
        "0.125\n",
        "phase_state_samples.csv");

    requireEqual(
        readAll(
            directory / "well_phase_rates.csv"),
        "phase_0_surface_rate_m3_s,phase_1_surface_rate_m3_s,phase_2_surface_rate_m3_s\n"
        "-2,-200,-6\n",
        "well_phase_rates.csv");

    requireEqual(
        readAll(
            directory / "adsorbed_mass_totals.csv"),
        "time_seconds,total_adsorbed_mass_kg,component_0_kg,component_1_kg,component_2_kg,component_3_kg,component_4_kg,component_5_kg\n"
        "12.5,21,1,2,3,4,5,6\n",
        "adsorbed_mass_totals.csv");

    requireEqual(
        readAll(
            directory / "trapped_mass_totals.csv"),
        "time_seconds,total_trapped_mass_kg,component_0_kg,component_1_kg,component_2_kg,component_3_kg,component_4_kg,component_5_kg\n"
        "12.5,21,1,2,3,4,5,6\n",
        "trapped_mass_totals.csv");

    requireEqual(
        readAll(
            directory / "trapped_mass_samples.csv"),
        "value_0,value_1\n"
        "1e-06,2e-06\n",
        "trapped_mass_samples.csv");

    std::filesystem::remove_all(directory);

    std::cout
        << "Output core validation: ALL PASS\n";

    return 0;
}
