/**
 * @file mesh_topology.cpp
 * @brief CanonicalMeshData 物化与 Mesh owner/neighbor 拓扑绑定。
 */
#include <cpgrid/mesh.hpp>

#include "mesh_canonical_data.hpp"

#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MPMC
{

void Mesh::initializeFromCanonical_(
    const detail::CanonicalMeshData &data)
{
    sourcePath_ = data.sourcePath;
    sourceFormat_ = data.sourceFormat;
    logicalDimensions_ = data.logicalDimensions;

    globalNodeCount_ = data.nodes.size();
    globalCellCount_ = data.cells.size();
    globalFaceInstanceCount_ = 0;
    for (const auto &cell : data.cells)
        globalFaceInstanceCount_ += cell.faces.size();
    globalUniqueFaceCount_ = data.faceNeighbors.size();

    nodes_.clear();
    cells_.clear();
    cartesianDirectory_.clear();
    faceNeighbors_.reserve(data.faceNeighbors.size());
    for (const auto &row : data.faceNeighbors)
    {
        faceNeighbors_.push_back(
            {{static_cast<PetscInt>(row[0]),
              static_cast<PetscInt>(row[1])}});
    }

    if (data.nodes.size() >
        static_cast<std::size_t>(
            std::numeric_limits<PetscInt>::max()))
    {
        throw std::overflow_error(
            "Mesh node count does not fit PetscInt.");
    }
    nodes_.reserve(data.nodes.size());
    for (std::size_t index = 0; index < data.nodes.size(); ++index)
    {
        const auto &node = data.nodes[index];
        nodes_.emplace_back(
            node.x,
            node.y,
            node.z,
            static_cast<PetscInt>(index));
    }

    cells_.reserve(data.cells.size());
    std::vector<PetscInt> cartesianByInput;
    cartesianByInput.reserve(data.cells.size());

    for (std::size_t storage = 0; storage < data.cells.size(); ++storage)
    {
        const auto &input = data.cells[storage];
        std::vector<Face> faces;
        faces.reserve(input.faces.size());

        for (const auto &inputFace : input.faces)
        {
            std::vector<Node *> faceNodes;
            faceNodes.reserve(inputFace.nodeIndices.size());
            for (std::size_t nodeIndex : inputFace.nodeIndices)
            {
                if (nodeIndex >= nodes_.size())
                    throw std::runtime_error(
                        "Canonical mesh face references invalid node id.");
                faceNodes.push_back(&nodes_[nodeIndex]);
            }

            faces.emplace_back(
                std::move(faceNodes),
                static_cast<PetscInt>(inputFace.inputPosition),
                static_cast<PetscInt>(inputFace.inputFaceId));
        }

        cells_.emplace_back(
            std::move(faces),
            static_cast<PetscInt>(storage));
        cartesianByInput.push_back(
            static_cast<PetscInt>(input.cartesianId));
    }

    cartesianDirectory_.reset(
        cartesianByInput,
        "Canonical mesh contains duplicate Cartesian cell id.");

    if (nodes_.empty() || cells_.empty())
        throw std::runtime_error(
            "Mesh data contain no nodes or cells.");
}

void Mesh::initializeTopology()
{
    if (topologyInitialized_)
        return;

    bindNeighbors_();
    topologyInitialized_ = true;
}

void Mesh::bindNeighbors_()
{
    std::vector<std::array<PetscInt, 2>> neighbors = faceNeighbors_;
    if (neighbors.empty())
    {
        const auto inputNeighbors =
            detail::loadMrstFaceNeighbors(dataDirectory_);
        neighbors.reserve(inputNeighbors.size());
        for (const auto &row : inputNeighbors)
        {
            neighbors.push_back(
                {{static_cast<PetscInt>(row[0]),
                  static_cast<PetscInt>(row[1])}});
        }
        faceNeighbors_ = neighbors;
    }
    globalUniqueFaceCount_ = neighbors.size();

    std::vector<std::array<Face *, 2>> faceInstances(neighbors.size());
    for (Polyhedron &cell : cells_)
    {
        cell.rebindFaceOwners();
        for (Face &face : cell.faces())
        {
            const PetscInt inputFaceId = face.inputFaceId();
            if (inputFaceId < 0 ||
                static_cast<std::size_t>(inputFaceId) >= neighbors.size())
            {
                throw std::runtime_error(
                    "G.faces.neighbors is inconsistent with face ids.");
            }

            auto &instances =
                faceInstances[static_cast<std::size_t>(inputFaceId)];
            if (instances[0] == nullptr)
                instances[0] = &face;
            else if (instances[1] == nullptr)
                instances[1] = &face;
            else
                throw std::runtime_error(
                    "A mesh face is referenced by more than two cells.");
        }
    }

    for (Polyhedron &cell : cells_)
    {
        const PetscInt storageId = cell.storageIndex();
        for (Face &face : cell.faces())
        {
            const PetscInt inputFaceId = face.inputFaceId();
            const auto &row =
                neighbors[static_cast<std::size_t>(inputFaceId)];
            const PetscInt first = row[0];
            const PetscInt second = row[1];
            PetscInt neighborStorage = -1;

            if (first == storageId)
                neighborStorage = second;
            else if (second == storageId)
                neighborStorage = first;
            else
                throw std::runtime_error(
                    "Face neighbor row does not contain its owner storage id.");

            if (neighborStorage < 0)
            {
                face.setNeighborCell(nullptr);
                face.setNeighborFace(nullptr);
                continue;
            }
            if (static_cast<std::size_t>(neighborStorage) >= cells_.size())
                throw std::runtime_error(
                    "Face references invalid neighbor storage id.");

            face.setNeighborCell(
                &cells_[static_cast<std::size_t>(neighborStorage)]);

            const auto &instances =
                faceInstances[static_cast<std::size_t>(inputFaceId)];
            Face *opposite =
                instances[0] == &face ? instances[1] : instances[0];
            if (opposite == nullptr ||
                opposite->ownerCell() != face.neighborCell())
            {
                throw std::runtime_error(
                    "An internal mesh face is missing its opposite cell-face instance.");
            }
            face.setNeighborFace(opposite);
        }
    }
}

} // namespace MPMC
