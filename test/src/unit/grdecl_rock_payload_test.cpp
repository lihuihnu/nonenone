/**
 * @file grdecl_rock_payload_test.cpp
 * @brief 验证 GRDECL expanded geometry 与 root-only 岩石载荷可以零拷贝分离。
 */
#include <cpgrid/grdecl_rock_payload.hpp>

#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

} // namespace

int main()
{
    MPMC::GrdeclGridData data;
    data.activeCells.resize(2);
    data.porosity = {0.1, 0.2};
    data.permeabilityX = {1.0, 2.0};
    data.permeabilityY = {3.0, 4.0};
    data.permeabilityZ = {5.0, 6.0};

    const auto payload = MPMC::takeRootRockData(data);

    require(payload.cellCount == 2, "active-cell count changed during rock extraction");
    require(payload.porosity == std::vector<double>({0.1, 0.2}), "porosity payload mismatch");
    require(payload.permeability[0] == std::vector<double>({1.0, 2.0}), "Kx payload mismatch");
    require(payload.permeability[1] == std::vector<double>({3.0, 4.0}), "Ky payload mismatch");
    require(payload.permeability[2] == std::vector<double>({5.0, 6.0}), "Kz payload mismatch");

    require(data.activeCells.size() == 2, "rock extraction must not mutate expanded geometry");
    require(data.porosity.empty(), "porosity source was not moved");
    require(data.permeabilityX.empty(), "Kx source was not moved");
    require(data.permeabilityY.empty(), "Ky source was not moved");
    require(data.permeabilityZ.empty(), "Kz source was not moved");

    std::cout << "GRDECL root rock payload extraction: ALL PASS\n";
    return 0;
}
