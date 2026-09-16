/**
 * @file cartesian_cell_directory_test.cpp
 * @brief 紧凑 Cartesian cell 目录的单调/乱序/重复与反向重建回归。
 */
#include <cpgrid/cartesian_cell_directory.hpp>

#include <cstdint>
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

template <class F>
void requireThrows(F &&callable, const std::string &expected)
{
    try
    {
        callable();
    }
    catch (const std::exception &caught)
    {
        require(caught.what() == expected, "unexpected exception text");
        return;
    }
    throw std::runtime_error("expected exception was not thrown");
}

} // namespace

int main()
{
    using Directory = MPMC::detail::CompactCartesianCellDirectory<std::int32_t>;

    Directory monotone;
    monotone.reset({2, 5, 9, 12}, "duplicate");
    require(monotone.size() == 4, "monotone size");
    require(monotone.usesImplicitInputOrder(), "monotone directory should be compact");
    require(monotone.inputIndex(2) == 0, "monotone first lookup");
    require(monotone.inputIndex(9) == 2, "monotone middle lookup");
    require(monotone.cartesianIdsInInputOrder() == std::vector<std::int32_t>({2, 5, 9, 12}),
            "monotone reverse reconstruction");

    Directory permuted;
    permuted.reset({9, 2, 12, 5}, "duplicate");
    require(!permuted.usesImplicitInputOrder(), "permuted directory needs a permutation");
    require(permuted.inputIndex(2) == 1, "permuted lookup 2");
    require(permuted.inputIndex(5) == 3, "permuted lookup 5");
    require(permuted.inputIndex(9) == 0, "permuted lookup 9");
    require(permuted.inputIndex(12) == 2, "permuted lookup 12");
    require(permuted.cartesianIdsInInputOrder() == std::vector<std::int32_t>({9, 2, 12, 5}),
            "permuted reverse reconstruction");

    requireThrows(
        [&] { (void)permuted.inputIndex(7); },
        "Cartesian cell id is not active in this Mesh.");

    Directory duplicate;
    requireThrows(
        [&] { duplicate.reset({1, 4, 4, 8}, "duplicate Cartesian"); },
        "duplicate Cartesian");

    std::cout << "[PASS] compact Cartesian cell directory\n";
    return 0;
}
