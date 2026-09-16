/**
 * @file grdecl_test.cpp
 * @brief 单元测试：验证 `grdecl` 的核心语义、边界条件和回归行为。
 */
#include <cpgrid/grdecl.hpp>
#include <cpgrid/csv_reader.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

bool close(double a, double b, double tolerance = 1.0e-12)
{
    return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

template <class Function>
void requireThrows(Function &&function, const std::string &message)
{
    try
    {
        function();
    }
    catch (const std::exception &)
    {
        return;
    }
    throw std::runtime_error(message);
}

void writeSingleCellDeck(
    const std::filesystem::path &path,
    const std::string &zcorn,
    const std::string &tail = {})
{
    std::ofstream deck(path);
    deck
        << "SPECGRID\n 1 1 1 1 F /\n"
        << "COORD\n"
        << " 0 0 0  0 0 10\n"
        << " 1 0 0  1 0 10\n"
        << " 0 1 0  0 1 10\n"
        << " 1 1 0  1 1 10 /\n"
        << "ZCORN\n " << zcorn << " /\n"
        << tail;
}

} // namespace

int main()
{
    try
    {
        const auto temp = std::filesystem::temp_directory_path() / "mpmc_grdecl_unit";
        std::filesystem::remove_all(temp);
        std::filesystem::create_directories(temp);

        {
            std::ofstream rock(temp / "rock.inc");
            rock << "PORO\n 0.20 0.30 /\n"
                 << "PERMX\n 100 200 /\n";
        }
        {
            std::ofstream deck(temp / "tiny.grdecl");
            deck
                << "-- direct GRDECL parser regression\n"
                << "METRIC\n"
                << "SPECGRID\n 2 1 1 1 F /\n"
                << "COORD\n"
                << " 0 0 0  0 0 10\n"
                << " 1 0 0  1 0 10\n"
                << " 2 0 0  2 0 10\n"
                << " 0 1 0  0 1 10\n"
                << " 1 1 0  1 1 10\n"
                << " 2 1 0  2 1 10 /\n"
                << "ZCORN\n 8*0 8*10 /\n"
                << "ACTNUM\n 2*1 /\n"
                << "INCLUDE\n 'rock.inc' /\n";
        }

        const auto data = MPMC::loadGrdecl(temp / "tiny.grdecl");
        require(data.nx == 2 && data.ny == 1 && data.nz == 1, "logical dimensions");
        require(data.cartesianCellCount == 2, "Cartesian count");
        require(data.activeCells.size() == 2, "active count");
        require(data.hasPorosity(), "PORO loaded");
        require(data.hasPermeability(), "PERM loaded");
        require(close(data.porosity[0], 0.20) && close(data.porosity[1], 0.30), "PORO values");

        constexpr double mdToM2 = 9.869232667160130e-16;
        require(std::abs(data.permeabilityX[0] - 100.0 * mdToM2) < 1.0e-28,
                "PERMX SI conversion");
        require(std::abs(data.permeabilityY[1] - 200.0 * mdToM2) < 1.0e-28,
                "PERMY fallback");
        require(std::abs(data.permeabilityZ[1] - 200.0 * mdToM2) < 1.0e-28,
                "PERMZ fallback");

        require(data.activeCells[0].neighborStorage[1] == 1, "I+ neighbor");
        require(data.activeCells[1].neighborStorage[0] == 0, "I- neighbor");
        require(data.activeCells[0].neighborStorage[0] == -1, "boundary neighbor");
        require(close(data.activeCells[0].corner[0].x, 0.0), "corner x0");
        require(close(data.activeCells[0].corner[7].x, 1.0), "corner x1");
        require(close(data.activeCells[0].corner[7].z, 10.0), "corner z bottom");

        MPMC::GrdeclLoadOptions raw;
        raw.convertToSi = false;
        raw.unitSystem = MPMC::GrdeclUnitSystem::Field;
        const auto fieldRaw = MPMC::loadGrdecl(temp / "tiny.grdecl", raw);
        require(!fieldRaw.convertedToSi, "raw conversion disabled");
        require(fieldRaw.unitSystem == MPMC::GrdeclUnitSystem::Field, "unit override");
        require(close(fieldRaw.permeabilityX[0], 100.0), "raw permeability preserved");

        MPMC::GrdeclLoadOptions fieldSi;
        fieldSi.unitSystem = MPMC::GrdeclUnitSystem::Field;
        const auto convertedField = MPMC::loadGrdecl(temp / "tiny.grdecl", fieldSi);
        require(close(convertedField.activeCells[0].corner[7].x, 0.3048),
                "FIELD geometry conversion");

        {
            std::ofstream deck(temp / "mapaxes.grdecl");
            deck
                << "GRIDUNIT\n METRES /\n"
                << "MAPUNITS\n METRES /\n"
                << "SPECGRID\n 1 1 1 1 F /\n"
                << "MAPAXES\n 5 5 10 10 15 5 /\n"
                << "COORD\n"
                << " 0 0 0  0 0 10\n"
                << " 5 0 0  5 0 10\n"
                << " 0 5 0  0 5 10\n"
                << " 5 5 0  5 5 10 /\n"
                << "ZCORN\n 4*0 4*10 /\n"
                << "PORO\n 0.20 /\n"
                << "PERMX\n 100 /\n";
        }
        const auto mapped = MPMC::loadGrdecl(temp / "mapaxes.grdecl");
        require(!mapped.unitSystemWasAssumed, "GRIDUNIT must define geometry units");
        require(close(mapped.activeCells[0].corner[0].x, 10.0) &&
                    close(mapped.activeCells[0].corner[0].y, 10.0),
                "MAPAXES origin translation");
        require(close(mapped.activeCells[0].corner[7].x, 10.0) &&
                    close(mapped.activeCells[0].corner[7].y, 10.0 - 5.0 * std::sqrt(2.0)),
                "MAPAXES rotation");

        {
            std::ofstream deck(temp / "inactive.grdecl");
            deck
                << "SPECGRID\n 2 1 1 1 F /\n"
                << "COORD\n"
                << " 0 0 0  0 0 10\n"
                << " 1 0 0  1 0 10\n"
                << " 2 0 0  2 0 10\n"
                << " 0 1 0  0 1 10\n"
                << " 1 1 0  1 1 10\n"
                << " 2 1 0  2 1 10 /\n"
                << "ZCORN\n 8*0 8*10 /\n"
                << "ACTNUM\n 1 0 /\n"
                << "PORO\n 0.20 0.30 /\n"
                << "PERMX\n 100 200 /\n";
        }
        const auto inactive = MPMC::loadGrdecl(temp / "inactive.grdecl");
        require(inactive.unitSystemWasAssumed, "missing unit system must be recorded");
        require(inactive.activeCells.size() == 1, "ACTNUM filtering");
        require(inactive.porosity.size() == 1 && close(inactive.porosity[0], 0.20),
                "active PORO remap");
        require(inactive.activeCells[0].neighborStorage[1] == -1,
                "inactive logical neighbor must become boundary");

        {
            std::ofstream deck(temp / "components.grdecl");
            deck
                << "SPECGRID\n 4 1 1 1 F /\n"
                << "COORD\n"
                << " 0 0 0  0 0 10\n"
                << " 1 0 0  1 0 10\n"
                << " 2 0 0  2 0 10\n"
                << " 3 0 0  3 0 10\n"
                << " 4 0 0  4 0 10\n"
                << " 0 1 0  0 1 10\n"
                << " 1 1 0  1 1 10\n"
                << " 2 1 0  2 1 10\n"
                << " 3 1 0  3 1 10\n"
                << " 4 1 0  4 1 10 /\n"
                << "ZCORN\n 16*0 16*10 /\n"
                << "ACTNUM\n 1 1 0 1 /\n"
                << "PORO\n 0.10 0.20 0.30 0.40 /\n"
                << "PERMX\n 10 20 30 40 /\n";
        }
        const auto allComponents = MPMC::loadGrdecl(temp / "components.grdecl");
        require(allComponents.activeCells.size() == 3,
                "default component policy must preserve all ACTNUM cells");
        MPMC::GrdeclLoadOptions largestOptions;
        largestOptions.keepLargestConnectedComponent = true;
        const auto largest = MPMC::loadGrdecl(temp / "components.grdecl", largestOptions);
        require(largest.activeCells.size() == 2, "largest component active count");
        require(largest.disconnectedCellCountRemoved == 1,
                "largest component removed count");
        require(largest.actnum == std::vector<int>({1, 1, 0, 0}),
                "largest component ACTNUM mask");
        require(largest.activeCells[0].cartesianId == 0 &&
                    largest.activeCells[1].cartesianId == 1,
                "largest component Cartesian ids");
        require(largest.porosity.size() == 2 && close(largest.porosity[1], 0.20),
                "largest component property remap");

        writeSingleCellDeck(
            temp / "trailing_numeric.grdecl",
            "4*0 4*10",
            "PORO\n 0.2invalid /\nPERMX\n 100 /\n");
        requireThrows(
            [&]() { (void)MPMC::loadGrdecl(temp / "trailing_numeric.grdecl"); },
            "GRDECL numeric tokens with trailing characters must be rejected");

        writeSingleCellDeck(
            temp / "oversized_repeat.grdecl",
            "999999999*0");
        requireThrows(
            [&]() { (void)MPMC::loadGrdecl(temp / "oversized_repeat.grdecl"); },
            "GRDECL repetition must be bounded before allocation");

        writeSingleCellDeck(
            temp / "non_finite.grdecl",
            "4*NaN 4*10");
        requireThrows(
            [&]() { (void)MPMC::loadGrdecl(temp / "non_finite.grdecl"); },
            "GRDECL non-finite geometry must be rejected");

        {
            std::ofstream cycle(temp / "cycle.grdecl");
            cycle << "INCLUDE\n 'cycle.grdecl' /\n";
        }
        requireThrows(
            [&]() { (void)MPMC::loadGrdecl(temp / "cycle.grdecl"); },
            "GRDECL INCLUDE cycles must be diagnosed immediately");

        {
            std::ofstream malformed(temp / "unterminated_quote.grdecl");
            malformed << "INCLUDE\n 'missing.inc /\n";
        }
        requireThrows(
            [&]() { (void)MPMC::loadGrdecl(temp / "unterminated_quote.grdecl"); },
            "GRDECL unterminated quotes must be rejected");

        {
            std::ofstream csv(temp / "invalid.csv");
            csv << "1invalid,2\n";
        }
        requireThrows(
            [&]() { (void)MPMC::readCsv<int>((temp / "invalid.csv").string()); },
            "CSV numeric tokens with trailing characters must be rejected");

        {
            std::ofstream csv(temp / "trailing_empty.csv");
            csv << "1,2,\n";
        }
        requireThrows(
            [&]() { (void)MPMC::readCsv<int>((temp / "trailing_empty.csv").string()); },
            "CSV trailing empty fields must be rejected");

        std::filesystem::remove_all(temp);
        std::cout << "GRDECL parser/corner-point topology: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "GRDECL TEST FAILED: " << error.what() << '\n';
        return 1;
    }
}
