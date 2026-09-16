/**
 * @file output_well_test.cpp
 * @brief 单元测试：验证 `output_well` 的核心语义、边界条件和回归行为。
 */
#include <output/output.hpp>

#include <indices/indices.hpp>
#include <well/well.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace
{

using Config =
    MPMC::CompositionalModelConfig<
        6, true, true,
        false, false, false>;

using Indices =
    MPMC::ScalarIndices<Config>;

using Well =
    MPMC::WellSpecification<
        Indices,
        int>;

std::string readAll(
    const std::filesystem::path &path)
{
    std::ifstream file(path);
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

} // namespace

int main()
{
    const auto directory =
        std::filesystem::temp_directory_path() /
        "mpmc_output_well_test";

    std::filesystem::remove_all(
        directory);

    MPMC::WellManager<Well>
        manager;

    manager.addWell(
        Well{
            10,
            "P1",
            MPMC::WellType::Producer,
            MPMC::WellControl::Bhp,
            1.0e7,
            0,
            {{0, 1.0}}});

    manager.addWell(
        Well{
            20,
            "I1",
            MPMC::WellType::Injector,
            MPMC::WellControl::Bhp,
            1.5e7,
            1,
            {{1, 1.0}}});

    std::vector<
        MPMC::WellState<Indices>>
        states(2);

    states[0].surfacePhaseRate =
        {-2.0, -200.0, -6.0};

    states[1].surfacePhaseRate =
        {0.0, 300.0, 0.0};

    {
        MPMC::OutputWriter<Indices>
            writer(directory);

        MPMC::writeWellStates(
            writer,
            manager,
            states,
            {20, 10});
    }

    const std::string expected =
        "phase_0_surface_rate_m3_s,phase_1_surface_rate_m3_s,phase_2_surface_rate_m3_s\n"
        "0,300,0\n"
        "-2,-200,-6\n";

    const std::string actual =
        readAll(
            directory / "well_phase_rates.csv");

    if (actual != expected)
    {
        throw std::runtime_error(
            "Well output ordering/format mismatch.");
    }

    std::filesystem::remove_all(
        directory);

    std::cout
        << "Output Well bridge: ALL PASS\n";

    return 0;
}
