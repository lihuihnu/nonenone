/**
 * @file lab_pseudo3d_geometry_test.cpp
 * @brief 回归：60x20x1 实验 slab 只锁定拓扑，尺寸由装置输入。
 */
#include "../../../case/scw_kerogen_common/lab_pseudo3d_slab_geometry.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}
}

int main()
{
    try
    {
        // Test-only dimensions: verify formulas, not a scientific apparatus choice.
        const ScwKerogenLab::Pseudo3DSlabGeometry g(0.30, 0.10, 0.01);
        require(g.nx == 60 && g.ny == 20 && g.nz == 1,
                "slab topology must remain 60x20x1");
        require(g.cellCount() == 1200,
                "slab cell count must be 1200");
        require(std::abs(g.dx() - 0.005) < 1.0e-15,
                "dx must equal Lx/60");
        require(std::abs(g.dy() - 0.005) < 1.0e-15,
                "dy must equal Ly/20");
        require(std::abs(g.dz() - 0.01) < 1.0e-15,
                "dz must equal physical slab thickness");
        require(std::abs(g.cellBulkVolume() -
                         g.bulkVolume() / 1200.0) < 1.0e-18,
                "uniform cell volume must close");

        bool rejected = false;
        try
        {
            (void)ScwKerogenLab::Pseudo3DSlabGeometry(0.30, 0.10, 0.0);
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
        }
        require(rejected, "unresolved/TBD zero slab thickness must be rejected");

        std::cout << "Lab pseudo3D slab geometry: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Lab pseudo3D slab geometry: FAIL: "
                  << error.what() << '\n';
        return 1;
    }
}
