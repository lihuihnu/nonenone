/**
 * @file grdecl_active.cpp
 * @brief GRDECL ACTNUM 连通域筛选与 active 属性压缩实现。
 */
#include "grdecl_detail.hpp"

#include <stdexcept>

namespace MPMC::grdecl_detail
{

std::size_t retainLargestConnectedComponent(
    std::vector<int> &actnum,
    int nx,
    int ny,
    int nz)
{
    std::vector<unsigned char> visited(actnum.size(), 0);
    std::vector<std::size_t> largest;
    std::vector<std::size_t> component;
    std::vector<std::size_t> pending;
    std::size_t activeCount = 0;

    for (std::size_t seed = 0; seed < actnum.size(); ++seed)
    {
        if (actnum[seed] == 0)
            continue;
        ++activeCount;
        if (visited[seed] != 0)
            continue;

        component.clear();
        pending.clear();
        pending.push_back(seed);
        visited[seed] = 1;
        for (std::size_t cursor = 0; cursor < pending.size(); ++cursor)
        {
            const std::size_t cart = pending[cursor];
            component.push_back(cart);
            const int i = static_cast<int>(cart % static_cast<std::size_t>(nx));
            const int j = static_cast<int>(
                (cart / static_cast<std::size_t>(nx)) % static_cast<std::size_t>(ny));
            const int k = static_cast<int>(
                cart / (static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny)));

            const auto visit = [&](int ni, int nj, int nk) {
                if (ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz)
                    return;
                const std::size_t neighbor = static_cast<std::size_t>(ni) +
                    static_cast<std::size_t>(nx) *
                        (static_cast<std::size_t>(nj) +
                         static_cast<std::size_t>(ny) * static_cast<std::size_t>(nk));
                if (actnum[neighbor] != 0 && visited[neighbor] == 0)
                {
                    visited[neighbor] = 1;
                    pending.push_back(neighbor);
                }
            };

            visit(i - 1, j, k);
            visit(i + 1, j, k);
            visit(i, j - 1, k);
            visit(i, j + 1, k);
            visit(i, j, k - 1);
            visit(i, j, k + 1);
        }

        // 数值：同规模连通分量按最小 Cartesian seed 确定，保证结果可复现。
        if (component.size() > largest.size())
            largest = component;
    }

    if (largest.empty() || largest.size() == activeCount)
        return 0;

    std::vector<unsigned char> keep(actnum.size(), 0);
    for (std::size_t cart : largest)
        keep[cart] = 1;
    for (std::size_t cart = 0; cart < actnum.size(); ++cart)
        if (actnum[cart] != 0 && keep[cart] == 0)
            actnum[cart] = 0;
    return activeCount - largest.size();
}

std::vector<double> activeProperty(
    const std::optional<std::vector<double>> &raw,
    const std::vector<int> &storageByCartesian,
    std::size_t cartesianCount,
    const std::string &keyword)
{
    if (!raw)
        return {};
    if (raw->size() != cartesianCount)
        throw std::runtime_error(keyword + " size does not match Cartesian cell count.");

    std::vector<double> out;
    out.reserve(cartesianCount);
    for (std::size_t cart = 0; cart < cartesianCount; ++cart)
        if (storageByCartesian[cart] >= 0)
            out.push_back((*raw)[cart]);
    return out;
}

} // namespace MPMC::grdecl_detail
