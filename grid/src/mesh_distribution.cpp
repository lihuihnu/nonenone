/**
 * @file mesh_distribution.cpp
 * @brief CpGrid root-only Mesh ingest 后的全局目录广播与 owned+ghost 拓扑分发。
 */
#include <cpgrid/distributed_mesh_loader.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace MPMC
{
namespace
{

constexpr int kSnapshotSizeTag = 17001;
constexpr int kSnapshotDataTag = 17002;
constexpr std::size_t kMaxMpiByteChunk =
    static_cast<std::size_t>(std::numeric_limits<int>::max());

class ByteWriter final
{
  public:
    template <class T>
    void pod(const T &value)
    {
        static_assert(std::is_trivially_copyable_v<T>);
        const auto *first = reinterpret_cast<const char *>(&value);
        data_.insert(data_.end(), first, first + sizeof(T));
    }

    void string(const std::string &value)
    {
        const auto size = static_cast<std::uint64_t>(value.size());
        pod(size);
        data_.insert(data_.end(), value.begin(), value.end());
    }

    template <class T>
    void vectorPod(const std::vector<T> &values)
    {
        static_assert(std::is_trivially_copyable_v<T>);
        const auto size = static_cast<std::uint64_t>(values.size());
        pod(size);
        if (values.empty())
            return;
        const auto *first = reinterpret_cast<const char *>(values.data());
        data_.insert(
            data_.end(),
            first,
            first + sizeof(T) * values.size());
    }

    [[nodiscard]] std::vector<char> release() &&
    {
        return std::move(data_);
    }

  private:
    std::vector<char> data_;
};

class ByteReader final
{
  public:
    explicit ByteReader(const std::vector<char> &data)
        : data_(data)
    {
    }

    template <class T>
    T pod()
    {
        static_assert(std::is_trivially_copyable_v<T>);
        require(sizeof(T));
        T value{};
        std::memcpy(&value, data_.data() + offset_, sizeof(T));
        offset_ += sizeof(T);
        return value;
    }

    std::string string()
    {
        const auto size = pod<std::uint64_t>();
        requireSize(size);
        const std::size_t count = static_cast<std::size_t>(size);
        require(count);
        std::string value(data_.data() + offset_, count);
        offset_ += count;
        return value;
    }

    template <class T>
    std::vector<T> vectorPod()
    {
        static_assert(std::is_trivially_copyable_v<T>);
        const auto size = pod<std::uint64_t>();
        requireSize(size);
        const std::size_t count = static_cast<std::size_t>(size);
        if (count > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::runtime_error("Distributed Mesh payload size overflows size_t.");
        const std::size_t bytes = count * sizeof(T);
        require(bytes);
        std::vector<T> values(count);
        if (bytes > 0)
            std::memcpy(values.data(), data_.data() + offset_, bytes);
        offset_ += bytes;
        return values;
    }

    void requireEnd() const
    {
        if (offset_ != data_.size())
            throw std::runtime_error("Distributed Mesh payload has trailing bytes.");
    }

  private:
    static void requireSize(std::uint64_t value)
    {
        if (value > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()))
            throw std::runtime_error("Distributed Mesh payload size exceeds size_t.");
    }

    void require(std::size_t bytes) const
    {
        if (bytes > data_.size() - offset_)
            throw std::runtime_error("Distributed Mesh payload is truncated.");
    }

    const std::vector<char> &data_;
    std::size_t offset_{0};
};

void broadcastBytes(
    CpCommunicator comm,
    int root,
    int rank,
    std::vector<char> &payload)
{
    unsigned long long size =
        rank == root ? static_cast<unsigned long long>(payload.size()) : 0ULL;
    PetscCallMPIAbort(
        comm,
        MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, root, comm));

    if (size > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()))
        throw std::runtime_error("Broadcast Mesh payload is too large for this process.");
    if (rank != root)
        payload.resize(static_cast<std::size_t>(size));

    std::size_t offset = 0;
    while (offset < payload.size())
    {
        const std::size_t remaining = payload.size() - offset;
        const int chunk = static_cast<int>(std::min(remaining, kMaxMpiByteChunk));
        PetscCallMPIAbort(
            comm,
            MPI_Bcast(payload.data() + offset, chunk, MPI_BYTE, root, comm));
        offset += static_cast<std::size_t>(chunk);
    }
}

void sendBytes(
    CpCommunicator comm,
    int target,
    const std::vector<char> &payload)
{
    const unsigned long long size =
        static_cast<unsigned long long>(payload.size());
    PetscCallMPIAbort(
        comm,
        MPI_Send(
            const_cast<unsigned long long *>(&size),
            1,
            MPI_UNSIGNED_LONG_LONG,
            target,
            kSnapshotSizeTag,
            comm));

    std::size_t offset = 0;
    while (offset < payload.size())
    {
        const std::size_t remaining = payload.size() - offset;
        const int chunk = static_cast<int>(std::min(remaining, kMaxMpiByteChunk));
        PetscCallMPIAbort(
            comm,
            MPI_Send(
                const_cast<char *>(payload.data() + offset),
                chunk,
                MPI_BYTE,
                target,
                kSnapshotDataTag,
                comm));
        offset += static_cast<std::size_t>(chunk);
    }
}

std::vector<char> receiveBytes(
    CpCommunicator comm,
    int root)
{
    unsigned long long size = 0ULL;
    MPI_Status status{};
    PetscCallMPIAbort(
        comm,
        MPI_Recv(
            &size,
            1,
            MPI_UNSIGNED_LONG_LONG,
            root,
            kSnapshotSizeTag,
            comm,
            &status));

    if (size > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()))
        throw std::runtime_error("Received Mesh snapshot is too large for this process.");
    std::vector<char> payload(static_cast<std::size_t>(size));

    std::size_t offset = 0;
    while (offset < payload.size())
    {
        const std::size_t remaining = payload.size() - offset;
        const int chunk = static_cast<int>(std::min(remaining, kMaxMpiByteChunk));
        PetscCallMPIAbort(
            comm,
            MPI_Recv(
                payload.data() + offset,
                chunk,
                MPI_BYTE,
                root,
                kSnapshotDataTag,
                comm,
                &status));
        offset += static_cast<std::size_t>(chunk);
    }
    return payload;
}

void broadcastRootStatus(
    CpCommunicator comm,
    int root,
    int rank,
    std::string &error)
{
    int success = rank == root && error.empty() ? 1 : 0;
    PetscCallMPIAbort(
        comm,
        MPI_Bcast(&success, 1, MPI_INT, root, comm));

    std::vector<char> errorBytes;
    if (rank == root)
        errorBytes.assign(error.begin(), error.end());
    broadcastBytes(comm, root, rank, errorBytes);

    if (success == 0)
    {
        if (rank != root)
            error.assign(errorBytes.begin(), errorBytes.end());
        throw std::runtime_error(
            "Root-only CpGrid ingest failed: " + error);
    }
}

std::int64_t toWire(PetscInt value) noexcept
{
    return static_cast<std::int64_t>(value);
}

PetscInt fromWire(std::int64_t value)
{
    if (value < static_cast<std::int64_t>(std::numeric_limits<PetscInt>::min()) ||
        value > static_cast<std::int64_t>(std::numeric_limits<PetscInt>::max()))
        throw std::overflow_error("Distributed Mesh id does not fit PetscInt.");
    return static_cast<PetscInt>(value);
}

} // namespace

std::unique_ptr<Mesh> Mesh::distributePreparedRoot(
    std::unique_ptr<Mesh> rootMesh,
    CpCommunicator comm,
    int root)
{
    int rank = 0;
    int size = 0;
    PetscCallMPIAbort(comm, MPI_Comm_rank(comm, &rank));
    PetscCallMPIAbort(comm, MPI_Comm_size(comm, &size));
    if (root < 0 || root >= size)
        throw std::out_of_range("Distributed Mesh root rank is outside communicator.");

    std::string preflightError;
    if (rank == root)
    {
        try
        {
            if (!rootMesh)
                throw std::invalid_argument("Root rank requires a prepared Mesh.");
            if (!rootMesh->isReadyForDofs())
                throw std::invalid_argument("Root Mesh must complete prepareForUse() before distribution.");
            if (rootMesh->rank() != root || rootMesh->processCount() != size)
                throw std::invalid_argument("Root Mesh communicator metadata do not match distribution communicator.");
            if (rootMesh->isLocalSnapshotCompacted())
                throw std::invalid_argument("Root Mesh must retain complete topology before distribution.");
        }
        catch (const std::exception &error)
        {
            preflightError = error.what();
        }
    }
    else if (rootMesh)
    {
        preflightError = "Non-root rank must not provide a root Mesh.";
    }

    /* Non-root local argument errors must also become collective. */
    int localArgumentOk = (rank == root || !rootMesh) ? 1 : 0;
    int allArgumentsOk = 0;
    PetscCallMPIAbort(
        comm,
        MPI_Allreduce(&localArgumentOk, &allArgumentsOk, 1, MPI_INT, MPI_MIN, comm));
    if (allArgumentsOk == 0 && rank == root && preflightError.empty())
        preflightError = "A non-root rank provided an unexpected root Mesh.";
    broadcastRootStatus(comm, root, rank, preflightError);

    std::unique_ptr<Mesh> result;
    if (rank == root)
        result = std::move(rootMesh);
    else
        result = std::unique_ptr<Mesh>(new Mesh(comm));

    /*
     * 全局只广播轻量编号目录。拓扑、节点与面只通过后续 point-to-point
     * snapshot 发送给实际需要它们的 rank。
     */
    std::vector<char> globalPayload;
    if (rank == root)
    {
        ByteWriter writer;
        writer.string(result->sourcePath_);
        writer.string(result->sourceFormat_);
        for (int value : result->logicalDimensions_)
            writer.pod<std::int32_t>(static_cast<std::int32_t>(value));
        writer.pod<std::uint64_t>(static_cast<std::uint64_t>(result->globalNodeCount_));
        writer.pod<std::uint64_t>(static_cast<std::uint64_t>(result->globalCellCount_));
        writer.pod<std::uint64_t>(static_cast<std::uint64_t>(result->globalFaceInstanceCount_));
        writer.pod<std::uint64_t>(static_cast<std::uint64_t>(result->globalUniqueFaceCount_));

        std::vector<std::int64_t> rankOffsets(result->rankCurrentOffsets_.size());
        std::transform(
            result->rankCurrentOffsets_.begin(),
            result->rankCurrentOffsets_.end(),
            rankOffsets.begin(),
            [](PetscInt value) { return toWire(value); });
        writer.vectorPod(rankOffsets);

        const auto cartesianByInput =
            result->cartesianDirectory_.cartesianIdsInInputOrder();
        std::vector<std::int64_t> cartesian(cartesianByInput.size());
        std::transform(
            cartesianByInput.begin(),
            cartesianByInput.end(),
            cartesian.begin(),
            [](PetscInt value) { return toWire(value); });
        writer.vectorPod(cartesian);

        std::vector<std::int64_t> current(result->currentIdByInputIndex_.size());
        std::transform(
            result->currentIdByInputIndex_.begin(),
            result->currentIdByInputIndex_.end(),
            current.begin(),
            [](PetscInt value) { return toWire(value); });
        writer.vectorPod(current);
        globalPayload = std::move(writer).release();
    }
    broadcastBytes(comm, root, rank, globalPayload);

    if (rank != root)
    {
        ByteReader reader(globalPayload);
        result->sourcePath_ = reader.string();
        result->sourceFormat_ = reader.string();
        for (int &value : result->logicalDimensions_)
            value = static_cast<int>(reader.pod<std::int32_t>());
        result->globalNodeCount_ = static_cast<std::size_t>(reader.pod<std::uint64_t>());
        result->globalCellCount_ = static_cast<std::size_t>(reader.pod<std::uint64_t>());
        result->globalFaceInstanceCount_ = static_cast<std::size_t>(reader.pod<std::uint64_t>());
        result->globalUniqueFaceCount_ = static_cast<std::size_t>(reader.pod<std::uint64_t>());

        const auto rankOffsets = reader.vectorPod<std::int64_t>();
        if (rankOffsets.size() != static_cast<std::size_t>(size) + 1)
            throw std::runtime_error("Distributed Mesh rank-offset directory size mismatch.");
        result->rankCurrentOffsets_.resize(rankOffsets.size());
        std::transform(
            rankOffsets.begin(),
            rankOffsets.end(),
            result->rankCurrentOffsets_.begin(),
            [](std::int64_t value) { return fromWire(value); });
        if (result->rankCurrentOffsets_.front() != 0 ||
            result->rankCurrentOffsets_.back() != static_cast<PetscInt>(result->globalCellCount_))
            throw std::runtime_error("Distributed Mesh rank-offset directory is inconsistent.");

        const auto cartesian = reader.vectorPod<std::int64_t>();
        const auto current = reader.vectorPod<std::int64_t>();
        reader.requireEnd();

        if (cartesian.size() != result->globalCellCount_ ||
            current.size() != result->globalCellCount_)
            throw std::runtime_error("Distributed Mesh global directory size mismatch.");

        std::vector<PetscInt> cartesianByInput(cartesian.size());
        result->currentIdByInputIndex_.resize(current.size());
        result->inputIndexByCurrentId_.assign(current.size(), PetscInt(-1));

        for (std::size_t input = 0; input < current.size(); ++input)
        {
            const PetscInt cartesianId = fromWire(cartesian[input]);
            const PetscInt currentId = fromWire(current[input]);
            cartesianByInput[input] = cartesianId;
            result->currentIdByInputIndex_[input] = currentId;
            if (currentId < 0 || static_cast<std::size_t>(currentId) >= current.size())
                throw std::runtime_error("Distributed Mesh current id is outside global directory.");
            auto &slot = result->inputIndexByCurrentId_[static_cast<std::size_t>(currentId)];
            if (slot >= 0)
                throw std::runtime_error("Distributed Mesh contains duplicate current id.");
            slot = static_cast<PetscInt>(input);
        }
        try
        {
            result->cartesianDirectory_.reset(
                cartesianByInput,
                "Distributed Mesh contains duplicate Cartesian id.");
        }
        catch (const std::logic_error &)
        {
            throw std::runtime_error(
                "Distributed Mesh contains duplicate Cartesian id.");
        }
    }

    std::vector<std::vector<PetscInt>> ownedByRank;
    if (rank == root)
    {
        ownedByRank.resize(static_cast<std::size_t>(size));
        for (const Polyhedron &cell : result->cells_)
        {
            const int owner = cell.processorId();
            if (owner < 0 || owner >= size)
                throw std::runtime_error(
                    "Distributed Mesh encountered an invalid owner rank.");
            ownedByRank[static_cast<std::size_t>(owner)].push_back(cell.id());
        }
        for (auto &owned : ownedByRank)
            std::sort(owned.begin(), owned.end());
    }

    auto makeSnapshot = [&](int target)
    {
        const auto &owned = ownedByRank.at(static_cast<std::size_t>(target));
        std::vector<PetscInt> ghost;
        for (PetscInt currentId : owned)
        {
            const Polyhedron &cell = result->cellByCurrentId(currentId);
            for (const Face &face : cell.faces())
            {
                const Polyhedron *neighbor = face.neighborCell();
                if (neighbor != nullptr && neighbor->processorId() != target)
                    ghost.push_back(neighbor->id());
            }
        }
        std::sort(ghost.begin(), ghost.end());
        ghost.erase(std::unique(ghost.begin(), ghost.end()), ghost.end());

        std::vector<PetscInt> selected = owned;
        selected.insert(selected.end(), ghost.begin(), ghost.end());

        std::vector<PetscInt> nodeIds;
        for (PetscInt currentId : selected)
        {
            const Polyhedron &cell = result->cellByCurrentId(currentId);
            for (const Face &face : cell.faces())
                for (const Node *node : face.nodes())
                {
                    if (node == nullptr)
                        throw std::logic_error("Distributed Mesh snapshot encountered a null node.");
                    nodeIds.push_back(node->id());
                }
        }
        std::sort(nodeIds.begin(), nodeIds.end());
        nodeIds.erase(std::unique(nodeIds.begin(), nodeIds.end()), nodeIds.end());

        ByteWriter writer;
        writer.pod<std::uint64_t>(static_cast<std::uint64_t>(owned.size()));
        for (PetscInt value : owned)
            writer.pod<std::int64_t>(toWire(value));
        writer.pod<std::uint64_t>(static_cast<std::uint64_t>(ghost.size()));
        for (PetscInt value : ghost)
            writer.pod<std::int64_t>(toWire(value));

        writer.pod<std::uint64_t>(static_cast<std::uint64_t>(nodeIds.size()));
        for (PetscInt nodeId : nodeIds)
        {
            const Node &node = result->nodes_.at(static_cast<std::size_t>(nodeId));
            writer.pod<std::int64_t>(toWire(node.id()));
            writer.pod<double>(node.x());
            writer.pod<double>(node.y());
            writer.pod<double>(node.z());
        }

        writer.pod<std::uint64_t>(static_cast<std::uint64_t>(selected.size()));
        for (PetscInt currentId : selected)
        {
            const Polyhedron &cell = result->cellByCurrentId(currentId);
            writer.pod<std::int64_t>(toWire(cell.id()));
            writer.pod<std::int64_t>(toWire(cell.inputIndex()));
            writer.pod<std::int32_t>(static_cast<std::int32_t>(cell.processorId()));
            writer.pod<std::uint64_t>(static_cast<std::uint64_t>(cell.faceCount()));
            for (const Face &face : cell.faces())
            {
                writer.pod<std::int64_t>(toWire(face.inputPosition()));
                writer.pod<std::int64_t>(toWire(face.inputFaceId()));
                const Polyhedron *neighbor = face.neighborCell();
                writer.pod<std::int64_t>(
                    neighbor == nullptr ? std::int64_t(-1) : toWire(neighbor->id()));
                writer.pod<std::uint64_t>(static_cast<std::uint64_t>(face.nodeCount()));
                for (const Node *node : face.nodes())
                    writer.pod<std::int64_t>(toWire(node->id()));
            }
        }
        return std::move(writer).release();
    };

    std::string distributionError;
    std::vector<char> snapshot;
    if (rank == root)
    {
        for (int target = 0; target < size; ++target)
        {
            if (target == root)
                continue;

            std::vector<char> targetPayload;
            if (distributionError.empty())
            {
                try
                {
                    targetPayload = makeSnapshot(target);
                }
                catch (const std::exception &caught)
                {
                    distributionError = caught.what();
                }
            }
            sendBytes(comm, target, targetPayload);
        }
    }
    else
    {
        snapshot = receiveBytes(comm, root);
    }

    /*
     * 所有 rank 都在反序列化前完成一次统一状态广播。即使 root 在构造某个
     * target snapshot 时发现拓扑错误，也会向所有 non-root 发出空消息并让
     * 全 communicator 一致抛错，避免部分 rank 已进入求解而其它 rank 挂起。
     */
    broadcastRootStatus(comm, root, rank, distributionError);
    if (rank == root)
        return result;

    ByteReader reader(snapshot);

    const auto ownedCount = reader.pod<std::uint64_t>();
    result->ownedCellIds_.resize(static_cast<std::size_t>(ownedCount));
    for (PetscInt &value : result->ownedCellIds_)
        value = fromWire(reader.pod<std::int64_t>());

    const auto ghostCount = reader.pod<std::uint64_t>();
    result->ghostCellIds_.resize(static_cast<std::size_t>(ghostCount));
    for (PetscInt &value : result->ghostCellIds_)
        value = fromWire(reader.pod<std::int64_t>());

    const auto nodeCount = reader.pod<std::uint64_t>();
    result->nodes_.clear();
    result->nodes_.reserve(static_cast<std::size_t>(nodeCount));
    std::unordered_map<PetscInt, std::size_t> nodeLocal;
    nodeLocal.reserve(static_cast<std::size_t>(nodeCount));
    for (std::uint64_t i = 0; i < nodeCount; ++i)
    {
        const PetscInt nodeId = fromWire(reader.pod<std::int64_t>());
        const double x = reader.pod<double>();
        const double y = reader.pod<double>();
        const double z = reader.pod<double>();
        nodeLocal.emplace(nodeId, result->nodes_.size());
        result->nodes_.emplace_back(x, y, z, nodeId);
    }

    struct NeighborRecord final
    {
        PetscInt currentId{-1};
        PetscInt faceId{-1};
    };

    const auto cellCount = reader.pod<std::uint64_t>();
    result->cells_.clear();
    result->cells_.reserve(static_cast<std::size_t>(cellCount));
    std::vector<std::vector<NeighborRecord>> neighbors;
    neighbors.reserve(static_cast<std::size_t>(cellCount));

    for (std::uint64_t c = 0; c < cellCount; ++c)
    {
        const PetscInt currentId = fromWire(reader.pod<std::int64_t>());
        const PetscInt inputIndex = fromWire(reader.pod<std::int64_t>());
        const int owner = static_cast<int>(reader.pod<std::int32_t>());
        const auto faceCount = reader.pod<std::uint64_t>();

        std::vector<Face> faces;
        std::vector<NeighborRecord> faceNeighbors;
        faces.reserve(static_cast<std::size_t>(faceCount));
        faceNeighbors.reserve(static_cast<std::size_t>(faceCount));
        for (std::uint64_t f = 0; f < faceCount; ++f)
        {
            const PetscInt inputPosition = fromWire(reader.pod<std::int64_t>());
            const PetscInt inputFaceId = fromWire(reader.pod<std::int64_t>());
            const PetscInt neighborId = fromWire(reader.pod<std::int64_t>());
            const auto faceNodeCount = reader.pod<std::uint64_t>();
            std::vector<Node *> faceNodes;
            faceNodes.reserve(static_cast<std::size_t>(faceNodeCount));
            for (std::uint64_t n = 0; n < faceNodeCount; ++n)
            {
                const PetscInt nodeId = fromWire(reader.pod<std::int64_t>());
                const auto found = nodeLocal.find(nodeId);
                if (found == nodeLocal.end())
                    throw std::runtime_error("Distributed Mesh face references a missing local node.");
                faceNodes.push_back(&result->nodes_[found->second]);
            }
            faces.emplace_back(std::move(faceNodes), inputPosition, inputFaceId);
            faceNeighbors.push_back({neighborId, inputFaceId});
        }

        result->cells_.emplace_back(std::move(faces), inputIndex);
        result->cells_.back().setId(currentId);
        result->cells_.back().setProcessorId(owner);
        neighbors.push_back(std::move(faceNeighbors));
    }
    reader.requireEnd();

    std::unordered_map<PetscInt, std::size_t> cellLocal;
    cellLocal.reserve(result->cells_.size());
    for (std::size_t local = 0; local < result->cells_.size(); ++local)
        cellLocal.emplace(result->cells_[local].id(), local);

    for (std::size_t local = 0; local < result->cells_.size(); ++local)
    {
        Polyhedron &cell = result->cells_[local];
        cell.rebindFaceOwners();
        for (std::size_t f = 0; f < cell.faces().size(); ++f)
        {
            Face &face = cell.faces()[f];
            const PetscInt neighborId = neighbors[local][f].currentId;
            if (neighborId < 0)
                continue;
            const auto found = cellLocal.find(neighborId);
            if (found == cellLocal.end())
            {
                if (cell.isOwnedBy(rank))
                    throw std::runtime_error("Distributed owned cell is missing a one-ring ghost neighbor.");
                continue;
            }
            Polyhedron &neighbor = result->cells_[found->second];
            face.setNeighborCell(&neighbor);
            Face *opposite = nullptr;
            for (Face &candidate : neighbor.faces())
                if (candidate.inputFaceId() == face.inputFaceId())
                {
                    opposite = &candidate;
                    break;
                }
            if (opposite == nullptr)
                throw std::runtime_error("Distributed Mesh cannot find opposite face instance.");
            face.setNeighborFace(opposite);
        }
    }

    result->localCellCount_ =
        result->rankCurrentOffsets_.at(static_cast<std::size_t>(rank) + 1) -
        result->rankCurrentOffsets_.at(static_cast<std::size_t>(rank));
    if (result->localCellCount_ != static_cast<PetscInt>(result->ownedCellIds_.size()))
        throw std::runtime_error("Distributed Mesh owned-cell count disagrees with rank offsets.");
    result->faceNeighbors_.clear();
    result->topologyInitialized_ = true;
    result->partitioned_ = true;
    result->currentIdsReady_ = true;
    result->localSnapshotCompacted_ = true;
    result->validateMaterializedOrdering_();
    return result;
}

DistributedMeshLoadResult loadDistributedMeshFromRoot(
    const std::string &inputPath,
    bool forceGrdecl,
    const GrdeclLoadOptions &grdeclOptions,
    CpCommunicator comm,
    int root)
{
    int rank = 0;
    int size = 0;
    PetscCallMPIAbort(comm, MPI_Comm_rank(comm, &rank));
    PetscCallMPIAbort(comm, MPI_Comm_size(comm, &size));
    if (root < 0 || root >= size)
        throw std::out_of_range("Distributed Mesh root rank is outside communicator.");

    std::unique_ptr<Mesh> rootMesh;
    std::optional<GrdeclGridData> rootGrdecl;
    DistributedMeshInputFormat format = DistributedMeshInputFormat::MrstCsv;
    std::string error;

    if (rank == root)
    {
        try
        {
            const std::filesystem::path input(inputPath);
            std::string extension = input.extension().string();
            std::transform(
                extension.begin(),
                extension.end(),
                extension.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            const bool directGrdecl =
                forceGrdecl ||
                (std::filesystem::is_regular_file(input) && extension == ".grdecl");

            if (directGrdecl)
            {
                format = DistributedMeshInputFormat::Grdecl;
                rootGrdecl = loadGrdecl(input, grdeclOptions);
                rootMesh = std::make_unique<Mesh>(*rootGrdecl, comm);
            }
            else
            {
                format = DistributedMeshInputFormat::MrstCsv;
                rootMesh = std::make_unique<Mesh>(inputPath, comm);
            }
            rootMesh->prepareForRootDistribution(root);
        }
        catch (const std::exception &caught)
        {
            error = caught.what();
        }
    }

    broadcastRootStatus(comm, root, rank, error);

    int formatValue = rank == root && format == DistributedMeshInputFormat::Grdecl ? 1 : 0;
    PetscCallMPIAbort(comm, MPI_Bcast(&formatValue, 1, MPI_INT, root, comm));
    format = formatValue == 1 ? DistributedMeshInputFormat::Grdecl : DistributedMeshInputFormat::MrstCsv;

    DistributedMeshLoadResult result;
    result.format = format;
    result.mesh = Mesh::distributePreparedRoot(std::move(rootMesh), comm, root);
    if (rank == root)
        result.grdeclOnRoot = std::move(rootGrdecl);
    return result;
}

} // namespace MPMC
