/**
 * @file grdecl_builder.cpp
 * @brief 由标准化 COORD/ZCORN/ACTNUM 构造 active corner-point topology。
 */
#include "grdecl_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace MPMC::grdecl_detail
{
namespace
{

std::size_t zcornIndex(
    int nx, int ny, int i, int j, int k,
    int ix, int iy, int iz)
{
    const std::size_t nx2 = static_cast<std::size_t>(2 * nx);
    const std::size_t ny2 = static_cast<std::size_t>(2 * ny);
    return static_cast<std::size_t>(2 * i + ix) +
           nx2 * (static_cast<std::size_t>(2 * j + iy) +
                  ny2 * static_cast<std::size_t>(2 * k + iz));
}

std::size_t pillarIndex(int nx, int i, int j)
{
    return 6 * (static_cast<std::size_t>(i) +
                static_cast<std::size_t>(nx + 1) * static_cast<std::size_t>(j));
}

GrdeclPoint cornerPoint(
    const std::vector<double> &coord,
    const std::vector<double> &zcorn,
    int nx,
    int ny,
    int i,
    int j,
    int k,
    int ix,
    int iy,
    int iz)
{
    const std::size_t p = pillarIndex(nx, i + ix, j + iy);
    const double xt = coord[p + 0];
    const double yt = coord[p + 1];
    const double zt = coord[p + 2];
    const double xb = coord[p + 3];
    const double yb = coord[p + 4];
    const double zb = coord[p + 5];
    const double z = zcorn[zcornIndex(nx, ny, i, j, k, ix, iy, iz)];

    double t = 0.0;
    const double dz = zb - zt;
    if (std::abs(dz) > 1.0e-14 * std::max({1.0, std::abs(zt), std::abs(zb)}))
        t = (z - zt) / dz;

    return GrdeclPoint{xt + t * (xb - xt), yt + t * (yb - yt), z};
}

} // namespace

ActiveGridBuild buildActiveGrid(
    const std::vector<double> &coord,
    const std::vector<double> &zcorn,
    const std::vector<int> &actnum,
    int nx,
    int ny,
    int nz)
{
    const std::size_t cartesianCount =
        static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz);

    ActiveGridBuild result;
    result.storageByCartesian.assign(cartesianCount, -1);
    result.activeCells.reserve(cartesianCount);

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const int cart = i + nx * (j + ny * k);
                if (actnum[static_cast<std::size_t>(cart)] == 0)
                    continue;

                GrdeclCell cell;
                cell.cartesianId = cart;
                cell.i = i;
                cell.j = j;
                cell.k = k;
                for (int iz = 0; iz < 2; ++iz)
                    for (int iy = 0; iy < 2; ++iy)
                        for (int ix = 0; ix < 2; ++ix)
                        {
                            const std::size_t corner =
                                static_cast<std::size_t>(ix + 2 * iy + 4 * iz);
                            cell.corner[corner] = cornerPoint(
                                coord, zcorn, nx, ny, i, j, k, ix, iy, iz);
                        }

                result.storageByCartesian[static_cast<std::size_t>(cart)] =
                    static_cast<int>(result.activeCells.size());
                result.activeCells.push_back(cell);
            }

    if (result.activeCells.empty())
        throw std::runtime_error("GRDECL contains no active cells after ACTNUM filtering.");

    const std::array<std::array<int, 3>, 6> direction{{
        {{-1, 0, 0}}, {{1, 0, 0}}, {{0, -1, 0}},
        {{0, 1, 0}}, {{0, 0, -1}}, {{0, 0, 1}}
    }};
    for (auto &cell : result.activeCells)
        for (std::size_t d = 0; d < direction.size(); ++d)
        {
            const int ni = cell.i + direction[d][0];
            const int nj = cell.j + direction[d][1];
            const int nk = cell.k + direction[d][2];
            if (ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz)
                continue;
            const int ncart = ni + nx * (nj + ny * nk);
            cell.neighborStorage[d] =
                result.storageByCartesian[static_cast<std::size_t>(ncart)];
        }

    return result;
}

} // namespace MPMC::grdecl_detail
