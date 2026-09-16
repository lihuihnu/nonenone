/**
 * @file metis_graph.cpp
 * @brief 由网格邻接关系构造 METIS 分区图的实现。
 */
#include <cpgrid/metis_graph.hpp>

#include <cpgrid/face.hpp>
#include <cpgrid/polyhedron.hpp>
#include <cpgrid/mesh.hpp>

#include <stdexcept>

namespace MPMC
{

METISCsrGraph
createMetisCsrGraph(
    const Mesh &mesh)
{
    if (!mesh.isTopologyInitialized())
    {
        throw std::logic_error(
            "METIS graph requires initialized Mesh topology.");
    }

    METISCsrGraph graph;

    graph.offsets.reserve(
        mesh.cellCount() + 1);
    graph.neighbors.reserve(mesh.faceInstanceCount());
    graph.offsets.push_back(0);

    for (std::size_t storage = 0;
         storage < mesh.cellCount();
         ++storage)
    {
        const Polyhedron &cell =
            mesh.cellByStorageIndex(storage);

        for (const Face &face : cell.faces())
        {
            const Polyhedron *neighbor =
                face.neighborCell();

            if (neighbor == nullptr)
            {
                continue;
            }

            graph.neighbors.push_back(
                static_cast<idx_t>(
                    neighbor->storageIndex()));
        }

        graph.offsets.push_back(
            static_cast<idx_t>(
                graph.neighbors.size()));
    }

    return graph;
}

} // namespace MPMC
