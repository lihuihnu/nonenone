/**
 * @file cpgrid_geometry_benchmark.cpp
 * @brief 对比 CpGrid 历史独立几何查询与融合 setup 几何内核的 CPU 成本。
 */
#include <cpgrid/face.hpp>
#include <cpgrid/node.hpp>
#include <cpgrid/polyhedron.hpp>

#include <array>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

namespace
{
struct CellFixture final
{
    CellFixture(PetscInt id, double x, double y, double z, double sx, double sy, double sz)
        : nodes{{
              MPMC::Node(x, y, z, id * 8 + 0),
              MPMC::Node(x + sx, y, z, id * 8 + 1),
              MPMC::Node(x + sx, y + sy, z, id * 8 + 2),
              MPMC::Node(x, y + sy, z, id * 8 + 3),
              MPMC::Node(x, y, z + sz, id * 8 + 4),
              MPMC::Node(x + sx, y, z + sz, id * 8 + 5),
              MPMC::Node(x + sx, y + sy, z + sz, id * 8 + 6),
              MPMC::Node(x, y + sy, z + sz, id * 8 + 7)}}
    {
        const auto node = [this](std::size_t index) { return &nodes[index]; };
        std::vector<MPMC::Face> faces;
        faces.reserve(6);
        faces.emplace_back(std::vector<MPMC::Node *>{node(0), node(3), node(2), node(1)}, 0, id * 6 + 0);
        faces.emplace_back(std::vector<MPMC::Node *>{node(4), node(5), node(6), node(7)}, 1, id * 6 + 1);
        faces.emplace_back(std::vector<MPMC::Node *>{node(0), node(1), node(5), node(4)}, 2, id * 6 + 2);
        faces.emplace_back(std::vector<MPMC::Node *>{node(3), node(7), node(6), node(2)}, 3, id * 6 + 3);
        faces.emplace_back(std::vector<MPMC::Node *>{node(0), node(4), node(7), node(3)}, 4, id * 6 + 4);
        faces.emplace_back(std::vector<MPMC::Node *>{node(1), node(2), node(6), node(5)}, 5, id * 6 + 5);
        cell = MPMC::Polyhedron(std::move(faces), id);
    }

    std::array<MPMC::Node, 8> nodes;
    MPMC::Polyhedron cell;
};

template <class Function>
double measure(Function &&function)
{
    const auto begin = std::chrono::steady_clock::now();
    function();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(end - begin).count();
}
} // namespace

int main()
{
    constexpr PetscInt count = 50000;
    std::vector<std::unique_ptr<CellFixture>> cells;
    cells.reserve(static_cast<std::size_t>(count));

    for (PetscInt id = 0; id < count; ++id)
    {
        const double sx = 1.0 + 1.0e-6 * static_cast<double>(id % 97);
        const double sy = 1.2 + 2.0e-6 * static_cast<double>(id % 83);
        const double sz = 0.8 + 3.0e-6 * static_cast<double>(id % 71);
        cells.push_back(std::make_unique<CellFixture>(
            id, 1.1 * static_cast<double>(id), 0.03 * static_cast<double>(id % 17),
            0.02 * static_cast<double>(id % 23), sx, sy, sz));
    }

    volatile double sink = 0.0;
    const double separate = measure([&]() {
        for (const auto &entry : cells)
        {
            const auto &cell = entry->cell;
            sink += cell.computeVolume();
            const MPMC::Point cellCentroid = cell.centroid();
            sink += cellCentroid.x() + cellCentroid.y() + cellCentroid.z();
            for (const MPMC::Face &face : cell.faces())
            {
                sink += face.computeArea();
                const MPMC::Point centroid = face.centroid();
                const MPMC::Point normal = face.unitNormal();
                sink += centroid.x() + centroid.y() + centroid.z();
                sink += normal.x() + normal.y() + normal.z();
            }
        }
    });

    const double fused = measure([&]() {
        for (const auto &entry : cells)
        {
            const auto &cell = entry->cell;
            const MPMC::PolyhedronGeometry cellGeometry = cell.computeGeometry();
            sink += cellGeometry.volume;
            sink += cellGeometry.centroid.x() + cellGeometry.centroid.y() + cellGeometry.centroid.z();
            for (const MPMC::Face &face : cell.faces())
            {
                const MPMC::PolygonGeometry geometry = face.computeGeometry();
                sink += geometry.area;
                sink += geometry.centroid.x() + geometry.centroid.y() + geometry.centroid.z();
                sink += geometry.unitNormal.x() + geometry.unitNormal.y() + geometry.unitNormal.z();
            }
        }
    });

    std::cout << "CpGrid geometry cells=" << count
              << " separate_s=" << separate
              << " fused_s=" << fused
              << " speedup=" << (separate / fused)
              << " sink=" << sink << '\n';
    return 0;
}
