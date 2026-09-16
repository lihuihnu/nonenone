/**
 * @file rock_property_validation_test.cpp
 * @brief 验证 CpGrid 传统与 root-only rock loader 共享的属性校验/零值均值规则。
 */
#include <cpgrid/rock_property_validation.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

template <class Function>
void requireError(Function &&function, const std::string &expected)
{
    try
    {
        function();
    }
    catch (const std::exception &error)
    {
        require(error.what() == expected, "Rock validation error text changed.");
        return;
    }
    throw std::runtime_error("Expected rock validation failure was not raised.");
}

} // namespace

int main()
{
    const std::vector<double> porosity{0.0, 0.2, 0.4};
    const std::array<std::vector<double>, 3> permeability{{
        {0.0, 2.0, 4.0},
        {1.0, 0.0, 5.0},
        {2.0, 4.0, 0.0}}};

    MPMC::CpGridRockLoadOptions options;
    options.replaceZeroPorosityWithMean = true;
    options.replaceZeroPermeabilityWithMean = true;
    const auto means = MPMC::detail::validateRockProperties(
        porosity, permeability, 3, options);

    require(std::abs(means.porosity - 0.3) < 1.0e-15,
            "Porosity replacement mean changed.");
    require(std::abs(means.permeability[0] - 3.0) < 1.0e-15,
            "Kx replacement mean changed.");
    require(std::abs(means.permeability[1] - 3.0) < 1.0e-15,
            "Ky replacement mean changed.");
    require(std::abs(means.permeability[2] - 3.0) < 1.0e-15,
            "Kz replacement mean changed.");

    requireError(
        [&]
        {
            MPMC::detail::validateRockProperties(
                porosity, permeability, 2, options);
        },
        "Rock property count must match active Mesh cell count.");

    auto invalidPorosity = porosity;
    invalidPorosity[0] = 1.1;
    requireError(
        [&]
        {
            MPMC::detail::validateRockProperties(
                invalidPorosity, permeability, 3, options);
        },
        "Porosity must be finite and inside [0,1].");

    auto invalidPermeability = permeability;
    invalidPermeability[1][0] = -1.0;
    requireError(
        [&]
        {
            MPMC::detail::validateRockProperties(
                porosity, invalidPermeability, 3, options);
        },
        "Permeability must be finite and non-negative.");

    requireError(
        [&]
        {
            MPMC::detail::validateRockProperties(
                std::vector<double>(3, 0.0), permeability, 3, options);
        },
        "Cannot replace zero porosity: no positive sample exists.");

    auto zeroAxis = permeability;
    zeroAxis[2].assign(3, 0.0);
    requireError(
        [&]
        {
            MPMC::detail::validateRockProperties(
                porosity, zeroAxis, 3, options);
        },
        "Cannot replace zero permeability: an axis has no positive sample.");

    std::cout << "CpGrid rock-property validation: ALL PASS\n";
    return 0;
}
