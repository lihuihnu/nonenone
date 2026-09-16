/**
 * @file well_core_test.cpp
 * @brief 单元测试：验证 `well_core` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <well/well.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
using Config = MPMC::CompositionalModelConfig<6, true, true>;
using Indices = MPMC::ScalarIndices<Config>;
using Well = MPMC::WellSpecification<Indices, int>;

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double value, double expected, double tolerance, const char *message)
{
    if (std::abs(value - expected) > tolerance)
        throw std::runtime_error(message);
}

void testSpecificationAndSchedule()
{
    Well producer(3, "P1", MPMC::WellType::Producer,
                  MPMC::WellControl::TotalRate, 12.0, 2,
                  {{2, 1.0e-12}, {7, 2.0e-12}});
    producer.schedule.openTime = 10.0;
    producer.schedule.closeTime = 20.0;
    producer.validate(10);

    near(producer.signedTarget(), -12.0, 0.0,
         "producer signed rate target");
    require(!producer.isActive(9.0), "well must be closed before open time");
    require(producer.isActive(10.0), "open time is inclusive");
    require(producer.isActive(20.0), "close time is inclusive");
    require(!producer.isActive(21.0), "well must be closed after close time");

    bool duplicateRejected = false;
    try
    {
        auto invalid = producer;
        invalid.perforations.push_back({2, 3.0e-12});
        invalid.validate(10);
    }
    catch (const std::invalid_argument &)
    {
        duplicateRejected = true;
    }
    require(duplicateRejected, "duplicate perforations must be rejected");

    Well injector(4, "I1", MPMC::WellType::Injector,
                  MPMC::WellControl::Bhp, 15.0e6, 1,
                  {{1, 1.0e-12}});
    injector.injectionPhaseFraction[Indices::Phase::water] = 1.0;
    injector.injectionComponentMassFraction[1] = 1.0;
    injector.validate(10);
    near(injector.signedTarget(), 15.0e6, 0.0,
         "BHP target must not be sign-converted");
}

void testManager()
{
    MPMC::WellManager<Well> manager;
    manager.addWell(Well(1, "A", MPMC::WellType::Producer,
                         MPMC::WellControl::Bhp, 10.0e6, 0,
                         {{0, 1.0e-12}}));
    manager.addWell(Well(2, "B", MPMC::WellType::Producer,
                         MPMC::WellControl::Bhp, 11.0e6, 1,
                         {{1, 1.0e-12}}));
    require(manager.size() == 2, "manager size");
    require(manager.atId(2).name == "B", "manager id lookup");

    bool duplicateRejected = false;
    try
    {
        manager.addWell(Well(2, "duplicate", MPMC::WellType::Producer,
                             MPMC::WellControl::Bhp, 11.0e6, 2,
                             {{2, 1.0e-12}}));
    }
    catch (const std::invalid_argument &)
    {
        duplicateRejected = true;
    }
    require(duplicateRejected, "duplicate well ids must be rejected");
}

void testVerticalBuilder()
{
    const auto cellId = [](int i, int j, int k)
    {
        return i + 10 * (j + 10 * k);
    };

    const std::vector<double> wi{1.0, 2.0, 3.0};
    const auto perforations =
        MPMC::makeVerticalPerforations<int>(2, 3, 4, wi, cellId);

    require(perforations.size() == 3, "vertical perforation count");
    require(perforations[0].currentCellId == 432,
            "vertical first stable cell id");
    require(perforations[2].currentCellId == 632,
            "vertical last stable cell id");
    near(perforations[1].wellIndex, 2.0, 0.0,
         "vertical explicit WI mapping");
}

void testPeaceman()
{
    MPMC::VerticalPeacemanCell cell;
    cell.dx = 10.0;
    cell.dy = 10.0;
    cell.completionLength = 5.0;
    cell.kx = 1.0e-13;
    cell.ky = 1.0e-13;
    cell.wellRadius = 0.1;
    cell.skin = 0.0;

    near(MPMC::peacemanEquivalentRadius(cell),
         1.9798989873223334, 1.0e-14,
         "isotropic square Peaceman radius");
    near(MPMC::verticalPeacemanWellIndex(cell),
         1.0522374459131186e-12, 1.0e-25,
         "vertical Peaceman WI");

    // Anisotropic result is invariant under simultaneous x/y and kx/ky swap.
    MPMC::VerticalPeacemanCell a = cell;
    a.dx = 20.0;
    a.dy = 8.0;
    a.kx = 2.0e-13;
    a.ky = 5.0e-14;
    MPMC::VerticalPeacemanCell b = a;
    std::swap(b.dx, b.dy);
    std::swap(b.kx, b.ky);
    near(MPMC::verticalPeacemanWellIndex(a),
         MPMC::verticalPeacemanWellIndex(b), 1.0e-25,
         "Peaceman x/y symmetry");
}

} // namespace

int main()
{
    testSpecificationAndSchedule();
    testManager();
    testVerticalBuilder();
    testPeaceman();
    std::cout << "Well core validation: ALL PASS\n";
    return 0;
}
