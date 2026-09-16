/**
 * @file cpgrid_geometry_kernel_test.cpp
 * @brief 验证 CpGrid 融合几何内核与历史独立几何查询 bit-level 等价。
 */
#include <cpgrid/face.hpp>
#include <cpgrid/node.hpp>
#include <cpgrid/polyhedron.hpp>

#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>

namespace
{
std::uint64_t bits(double value)
{
    std::uint64_t result = 0;
    static_assert(sizeof(result) == sizeof(value));
    std::memcpy(&result, &value, sizeof(value));
    return result;
}

bool samePoint(const MPMC::Point &lhs, const MPMC::Point &rhs)
{
    return bits(lhs.x()) == bits(rhs.x()) &&
           bits(lhs.y()) == bits(rhs.y()) &&
           bits(lhs.z()) == bits(rhs.z());
}

struct CellFixture final
{
    explicit CellFixture(PetscInt id, double offset)
        : nodes{{
              MPMC::Node(0.00 + offset, 0.00, 0.00, id * 8 + 0),
              MPMC::Node(1.10 + offset, 0.03, 0.01, id * 8 + 1),
              MPMC::Node(1.00 + offset, 1.20, 0.07, id * 8 + 2),
              MPMC::Node(-0.06 + offset, 1.00, -0.02, id * 8 + 3),
              MPMC::Node(0.04 + offset, -0.02, 0.90, id * 8 + 4),
              MPMC::Node(1.13 + offset, 0.04, 1.02, id * 8 + 5),
              MPMC::Node(1.02 + offset, 1.18, 1.11, id * 8 + 6),
              MPMC::Node(-0.03 + offset, 1.04, 0.97, id * 8 + 7)}}
    {
        const auto node = [this](std::size_t index)
        {
            return &nodes[index];
        };

        std::vector<MPMC::Face> faces;
        faces.reserve(6);
        faces.emplace_back(
            std::vector<MPMC::Node *>{node(0), node(3), node(2), node(1)}, 0, id * 6 + 0);
        faces.emplace_back(
            std::vector<MPMC::Node *>{node(4), node(5), node(6), node(7)}, 1, id * 6 + 1);
        faces.emplace_back(
            std::vector<MPMC::Node *>{node(0), node(1), node(5), node(4)}, 2, id * 6 + 2);
        faces.emplace_back(
            std::vector<MPMC::Node *>{node(3), node(7), node(6), node(2)}, 3, id * 6 + 3);
        faces.emplace_back(
            std::vector<MPMC::Node *>{node(0), node(4), node(7), node(3)}, 4, id * 6 + 4);
        faces.emplace_back(
            std::vector<MPMC::Node *>{node(1), node(2), node(6), node(5)}, 5, id * 6 + 5);

        cell = MPMC::Polyhedron(std::move(faces), id);
    }

    std::array<MPMC::Node, 8> nodes;
    MPMC::Polyhedron cell;
};
} // namespace

int main()
{
    for (PetscInt id = 0; id < 1000; ++id)
    {
        CellFixture fixture(id, static_cast<double>(id) * 1.0e-5);

        const MPMC::PolyhedronGeometry fusedCell = fixture.cell.computeGeometry();
        if (bits(fusedCell.volume) != bits(fixture.cell.computeVolume()) ||
            !samePoint(fusedCell.centroid, fixture.cell.centroid()))
        {
            std::cerr << "CpGrid fused cell geometry changed historical arithmetic.\n";
            return 1;
        }

        for (const MPMC::Face &face : fixture.cell.faces())
        {
            const MPMC::PolygonGeometry fusedFace = face.computeGeometry();
            if (bits(fusedFace.area) != bits(face.computeArea()) ||
                !samePoint(fusedFace.centroid, face.centroid()) ||
                !samePoint(fusedFace.unitNormal, face.unitNormal()))
            {
                std::cerr << "CpGrid fused face geometry changed historical arithmetic.\n";
                return 2;
            }
        }
    }

    std::cout << "CpGrid fused geometry kernel: PASS\n";
    return 0;
}
