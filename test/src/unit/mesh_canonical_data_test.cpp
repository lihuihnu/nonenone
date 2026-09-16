/**
 * @file mesh_canonical_data_test.cpp
 * @brief 验证 MRST/GRDECL 到 CanonicalMeshData 的统一输入适配层。
 */
#include "mesh_canonical_data.hpp"

#include <cpgrid/grdecl.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void writeFile(
    const std::filesystem::path &filename,
    const char *content)
{
    std::filesystem::create_directories(filename.parent_path());
    std::ofstream file(filename);
    if (!file)
        throw std::runtime_error("Cannot create canonical-mesh test input.");
    file << content;
}

std::filesystem::path makeMrstFixture()
{
    const auto root =
        std::filesystem::temp_directory_path() /
        "mpmc_mesh_canonical_data_test";
    std::filesystem::remove_all(root);

    writeFile(root / "G/nodes/coords/data.csv",
              "0,0,0\n1,0,0\n0,1,0\n1,1,0\n"
              "0,0,1\n1,0,1\n0,1,1\n1,1,1\n");
    writeFile(root / "G/cells/indexMap/data.csv", "1\n");
    writeFile(root / "G/cells/faces/data.csv",
              "1,1\n2,2\n3,3\n4,4\n5,5\n6,6\n");
    writeFile(root / "G/cells/facePos/data.csv", "1\n7\n");
    writeFile(root / "G/faces/nodes/data.csv",
              "1\n3\n7\n5\n2\n6\n8\n4\n"
              "1\n5\n6\n2\n3\n4\n8\n7\n"
              "1\n2\n4\n3\n5\n7\n8\n6\n");
    writeFile(root / "G/faces/nodePos/data.csv",
              "1\n5\n9\n13\n17\n21\n25\n");
    writeFile(root / "G/faces/neighbors/data.csv",
              "1,0\n1,0\n1,0\n1,0\n1,0\n1,0\n");

    return root;
}

MPMC::GrdeclCell makeCell(
    int cartesianId,
    double x0,
    double x1)
{
    MPMC::GrdeclCell cell;
    cell.cartesianId = cartesianId;
    for (int iz = 0; iz < 2; ++iz)
        for (int iy = 0; iy < 2; ++iy)
            for (int ix = 0; ix < 2; ++ix)
            {
                const int corner = ix + 2 * iy + 4 * iz;
                cell.corner[static_cast<std::size_t>(corner)] = {
                    ix ? x1 : x0,
                    static_cast<double>(iy),
                    static_cast<double>(iz)};
            }
    return cell;
}

void testMrst()
{
    const auto root = makeMrstFixture();
    const auto data =
        MPMC::detail::loadMrstCanonicalMeshData(root.string());

    require(data.sourceFormat == "MRST_CSV", "MRST source format mismatch.");
    require(data.nodes.size() == 8, "MRST node count mismatch.");
    require(data.cells.size() == 1, "MRST cell count mismatch.");
    require(data.cells[0].cartesianId == 0, "MRST Cartesian id mismatch.");
    require(data.cells[0].faces.size() == 6, "MRST face count mismatch.");
    require(data.cells[0].faces[0].nodeIndices ==
                std::vector<std::size_t>({0, 2, 6, 4}),
            "MRST face-node ordering mismatch.");

    const auto neighbors =
        MPMC::detail::loadMrstFaceNeighbors(root.string());
    require(neighbors.size() == 6, "MRST neighbor count mismatch.");
    require(neighbors[0][0] == 0 && neighbors[0][1] == -1,
            "MRST boundary neighbor conversion mismatch.");

    std::filesystem::remove_all(root);
}

void testGrdecl()
{
    MPMC::GrdeclGridData input;
    input.sourcePath = "canonical_test.grdecl";
    input.nx = 2;
    input.ny = 1;
    input.nz = 1;
    input.nodeMergeTolerance = 1.0e-9;

    auto first = makeCell(0, 0.0, 1.0);
    auto second = makeCell(1, 1.0, 2.0);
    first.neighborStorage = {{-1, 1, -1, -1, -1, -1}};
    second.neighborStorage = {{0, -1, -1, -1, -1, -1}};
    input.activeCells = {first, second};

    const auto data =
        MPMC::detail::makeGrdeclCanonicalMeshData(input);

    require(data.sourceFormat == "GRDECL", "GRDECL source format mismatch.");
    require(data.logicalDimensions == std::array<int, 3>{{2, 1, 1}},
            "GRDECL dimensions mismatch.");
    require(data.nodes.size() == 12,
            "GRDECL conforming corners were not merged deterministically.");
    require(data.cells.size() == 2, "GRDECL cell count mismatch.");
    require(data.faceNeighbors.size() == 11, "GRDECL unique face count mismatch.");
    require(data.cells[0].faces[1].inputFaceId ==
                data.cells[1].faces[0].inputFaceId,
            "GRDECL shared face id mismatch.");
    require(data.cells[0].faces[1].nodeIndices ==
                std::vector<std::size_t>({1, 5, 7, 3}),
            "GRDECL face-node ordering mismatch.");
}

} // namespace

int main()
{
    testMrst();
    testGrdecl();
    std::cout << "[PASS] mesh canonical input adapters\n";
    return 0;
}
