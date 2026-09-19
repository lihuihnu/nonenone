#pragma once

#include <cmath>
#include <stdexcept>

/**
 * @file lab_pseudo3d_slab_geometry.hpp
 * @brief 实验装置驱动的 60x20x1 pseudo-3D slab 几何契约。
 *
 * 只锁定离散拓扑，不锁定任何物理尺寸。Lx/Ly/Lz 必须来自实际高温高压
 * slab/core-holder 的有效流动尺寸；旧 1.2m x 0.1m x 0.1m benchmark
 * 尺寸不得作为默认值进入新的实验可复现算例。
 */
namespace ScwKerogenLab
{

struct Pseudo3DSlabGeometry
{
    static constexpr int nx = 60;
    static constexpr int ny = 20;
    static constexpr int nz = 1;

    double lx{0.0};
    double ly{0.0};
    double lz{0.0};

    Pseudo3DSlabGeometry(double length, double inPlaneWidth, double thickness)
        : lx(length), ly(inPlaneWidth), lz(thickness)
    {
        if (!(lx > 0.0) || !(ly > 0.0) || !(lz > 0.0) ||
            !std::isfinite(lx) || !std::isfinite(ly) || !std::isfinite(lz))
        {
            throw std::invalid_argument(
                "Laboratory slab dimensions must be finite and positive.");
        }
    }

    [[nodiscard]] constexpr int cellCount() const noexcept
    {
        return nx * ny * nz;
    }

    [[nodiscard]] double dx() const noexcept { return lx / nx; }
    [[nodiscard]] double dy() const noexcept { return ly / ny; }
    [[nodiscard]] double dz() const noexcept { return lz; }

    [[nodiscard]] double bulkVolume() const noexcept
    {
        return lx * ly * lz;
    }

    [[nodiscard]] double cellBulkVolume() const noexcept
    {
        return bulkVolume() / static_cast<double>(cellCount());
    }
};

} // namespace ScwKerogenLab
