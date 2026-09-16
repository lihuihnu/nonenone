/**
 * @file peaceman.hpp
 * @brief Peaceman 井指数及等效井半径计算。
 */
#pragma once

#include <cmath>
#include <stdexcept>

namespace MPMC
{

struct VerticalPeacemanCell final
{
    double dx{0.0};
    double dy{0.0};
    double completionLength{0.0};
    double kx{0.0};
    double ky{0.0};
    double wellRadius{0.0};
    double skin{0.0};
};

/**
 * @brief 各向异性矩形笛卡尔单元中居中直井的 Peaceman 等效半径。
 *
 * 当前井模型采用以下 Peaceman 公式：
 *
 * r_e = 0.28 * sqrt(dx^2*sqrt(ky/kx) + dy^2*sqrt(kx/ky))
 *       / ((ky/kx)^(1/4) + (kx/ky)^(1/4)).
 */
[[nodiscard]] inline double peacemanEquivalentRadius(
    const VerticalPeacemanCell &cell)
{
    if (!(cell.dx > 0.0) || !(cell.dy > 0.0) ||
        !(cell.kx > 0.0) || !(cell.ky > 0.0) ||
        !std::isfinite(cell.dx) || !std::isfinite(cell.dy) ||
        !std::isfinite(cell.kx) || !std::isfinite(cell.ky))
    {
        throw std::invalid_argument(
            "Peaceman radius requires positive finite dx, dy, kx and ky.");
    }

    const double kyOverKx = cell.ky / cell.kx;
    const double kxOverKy = cell.kx / cell.ky;

    const double numerator =
        0.28 * std::sqrt(
                   cell.dx * cell.dx * std::sqrt(kyOverKx) +
                   cell.dy * cell.dy * std::sqrt(kxOverKy));

    const double denominator =
        std::pow(kyOverKx, 0.25) +
        std::pow(kxOverKy, 0.25);

    return numerator / denominator;
}

/**
 * @brief MPMC 直井模型使用的 Peaceman 产能指数。
 *
 * WI = 2*pi*h*sqrt(kx*ky) / (ln(re/rw) + skin).
 *
 * 渗透率必须与周围 Darcy 离散使用同一套内部 SI 单位。
 */
[[nodiscard]] inline double verticalPeacemanWellIndex(
    const VerticalPeacemanCell &cell)
{
    if (!(cell.completionLength > 0.0) ||
        !(cell.wellRadius > 0.0) ||
        !std::isfinite(cell.completionLength) ||
        !std::isfinite(cell.wellRadius) ||
        !std::isfinite(cell.skin))
    {
        throw std::invalid_argument(
            "Peaceman WI requires positive finite completion length/radius and finite skin.");
    }

    const double equivalentRadius =
        peacemanEquivalentRadius(cell);

    const double denominator =
        std::log(equivalentRadius / cell.wellRadius) + cell.skin;

    if (!(denominator > 0.0) || !std::isfinite(denominator))
        throw std::invalid_argument(
            "Peaceman WI denominator must be positive and finite; check grid size, well radius and skin.");

    constexpr double pi = 3.141592653589793238462643383279502884;
    const double effectivePermeability = std::sqrt(cell.kx * cell.ky);

    return 2.0 * pi * cell.completionLength * effectivePermeability /
           denominator;
}

} // namespace MPMC
