#ifndef DISKPACK_CORONA_H
#define DISKPACK_CORONA_H

#include <diskpack/geometry.h>
#include <diskpack/operator_table.h>
#include <vector>
#include <deque>
#include <memory>
#include <optional>

namespace diskpack {

class Corona;
class CoronaSignature;
class ConnectivityGraph;

using DiskPointer = std::shared_ptr<Disk>;
using CoronaSignaturePointer = std::shared_ptr<CoronaSignature>;

/// Stores information about disk adjacency patterns
class CoronaSignature {
public:
    explicit CoronaSignature(const Corona& corona);

    size_t baseType() const;
    bool testRadii(OperatorLookupTable& table) const;
    bool operator<(const CoronaSignature& other) const;
    bool operator==(const CoronaSignature& other) const;

private:
    const size_t numRadii;
    const size_t baseDiskType;
    std::vector<size_t> adjacencyCounts;
    std::vector<size_t> diskIndices;
    std::vector<bool> diskPresence;

    size_t& getAdjacency(size_t i, size_t j);
    size_t getAdjacency(size_t i, size_t j) const;
};

/// Manages disk connectivity patterns
class ConnectivityGraph {
public:
    static constexpr size_t MAX_SIGNATURES = 5000;

    explicit ConnectivityGraph(OperatorLookupTable& table);

    size_t size() const;
    void display() const;
    bool isViable() const;
    bool hasTriangle(size_t i, size_t j, size_t k) const;
    void refine(OperatorLookupTable& table);
    bool restore();

private:
    using SignatureList = std::list<CoronaSignaturePointer>;
    using VersionStack = std::stack<std::shared_ptr<SignatureList>>;

    BaseType precisionThreshold;
    bool overflowFlag = false;
    std::queue<std::tuple<size_t, size_t, size_t>> invalidTriangles;
    std::vector<SignatureList> signaturesByBase;
    std::vector<VersionStack> versionHistory;
    std::vector<std::vector<size_t>> transitionCounts;
    std::vector<std::vector<bool>> adjacencyMatrix;

    void addSignature(const CoronaSignature& signature);
    void removeSignature(const CoronaSignature& signature);
    void updateAdjacency();
    void processInvalidTriangles(
        std::shared_ptr<std::vector<std::shared_ptr<SignatureList>>> diff = nullptr
    );
    void fillCorona(Corona& corona, std::vector<DiskPointer>& packing,
                   size_t startIndex, std::set<CoronaSignaturePointer>& uniqueSignatures);

    size_t& getTransition(size_t base, size_t i, size_t j);
    size_t getTransition(size_t base, size_t i, size_t j) const;
};

/// Represents a ring of disks around central disk
class Corona {
public:
    Corona(const Disk& base, const std::vector<DiskPointer>& packing,
          OperatorLookupTable& table);

    bool isComplete() const;
    bool isContinuous() const;
    const Disk& baseDisk() const;

    struct Placement {
        bool success;
        Disk disk;

        Placement(bool s, const Disk& d) : success(s), disk(d) {}
    };

    Placement calculatePlacement(
        size_t diskType,
        const std::optional<ConnectivityGraph>& graph = std::nullopt
    );
    void addDisk(const DiskPointer& disk, size_t diskType);
    void removeDisk();
    void display() const;

private:
    static constexpr size_t DEFAULT_OPERATORS = 12;

    const Disk& centralDisk;
    OperatorLookupTable& opTable;
    std::deque<DiskPointer> ringDisks;
    std::vector<SSORef> frontOperators;
    std::vector<SSORef> backOperators;
    IntervalPair frontLeaf;
    IntervalPair backLeaf;
    std::stack<bool> operationHistory;

    SpiralSimilarityOperator composeOperators(
        size_t start, size_t end,
        const std::vector<SSORef>& ops
    ) const;
    void sortDisks();

    friend class CoronaSignature;
    friend class ConnectivityGraph;
};

} // namespace diskpack

#endif // DISKPACK_CORONA_H