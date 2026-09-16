/**
 * @file structured_cartesian_geometry_test.cpp
 * @brief StructuredGrid 可推导 Cartesian 几何量的 PETSc-free 回归测试。
 */
#include <structuredgrid/structured_cartesian_geometry.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main()
{
    using MPMC::StructuredFaceOrder;
    using namespace MPMC::detail;

    const std::array<double, 3> widths{2.0, 3.0, 5.0};
    assert(structuredCellVolume(widths) == 30.0);

    assert(structuredFaceHalfDistance(widths, StructuredFaceOrder::XMinus) == 1.0);
    assert(structuredFaceHalfDistance(widths, StructuredFaceOrder::YPlus) == 1.5);
    assert(structuredFaceHalfDistance(widths, StructuredFaceOrder::ZMinus) == 2.5);

    assert(structuredFaceNormalComponent(widths, StructuredFaceOrder::XMinus) == -15.0);
    assert(structuredFaceNormalComponent(widths, StructuredFaceOrder::XPlus) == 15.0);
    assert(structuredFaceNormalComponent(widths, StructuredFaceOrder::YMinus) == -10.0);
    assert(structuredFaceNormalComponent(widths, StructuredFaceOrder::YPlus) == 10.0);
    assert(structuredFaceNormalComponent(widths, StructuredFaceOrder::ZMinus) == -6.0);
    assert(structuredFaceNormalComponent(widths, StructuredFaceOrder::ZPlus) == 6.0);

    const auto centers = structuredAxisCenters(std::vector<double>{0.1, 0.2, 0.4});
    assert(centers.size() == 3);
    assert(centers[0] == 0.05);
    assert(centers[1] == 0.2);
    assert(centers[2] == 0.5);

    bool invalidDirection = false;
    try
    {
        (void)structuredFaceHalfDistance(widths, StructuredFaceOrder::Count);
    }
    catch (const std::invalid_argument &)
    {
        invalidDirection = true;
    }
    assert(invalidDirection);

    std::cout << "Structured Cartesian derived geometry: ALL PASS\n";
    return 0;
}
