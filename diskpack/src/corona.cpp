// diskpack/src/corona.cpp
#include "diskpack/corona.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <set>
#include <stdexcept>
#include <numeric>

namespace diskpack {

// -------------------------------
// Helpers
// -------------------------------

static IntervalPoint makeIntervalPoint(double x, double y, Interval uncert = Interval(1e-9)) {
    // create symmetric uncertainty around given coordinates
    double e = uncert.upper();
    return { Interval(x - e, x + e), Interval(y - e, y + e) };
}

// map triangular index -> pair (i,j)
// inverse of idx = (j*(j+1))/2 + i  with 0 <= i <= j
static inline void triangularIndexToPair(size_t idx, size_t &out_i, size_t &out_j) {
    // binary search j such that (j*(j+1))/2 <= idx < ((j+1)*(j+2))/2
    size_t j = static_cast<size_t>(std::floor((std::sqrt(8.0 * double(idx) + 1.0) - 1.0) / 2.0));
    while ((j * (j + 1)) / 2 > idx) --j;
    while (((j + 1) * (j + 2)) / 2 <= idx) ++j;
    size_t i = idx - (j * (j + 1)) / 2;
    out_i = i;
    out_j = j;
}

// -------------------------------
// CoronaSignature
// -------------------------------

CoronaSignature::CoronaSignature(const Corona& corona)
    : numRadii(corona.opTable.radii.size()),
      baseDiskType(corona.centralDisk.type()),
      adjacencyCounts((numRadii * numRadii + numRadii) / 2, 0),
      diskIndices(corona.ringDisks.size()),
      diskPresence(numRadii, false)
{
    if (!corona.isComplete()) {
        throw std::runtime_error("CoronaSignature: corona must be complete");
    }

    // extract disk types from the corona
    for (size_t k = 0; k < corona.ringDisks.size(); ++k) {
        diskIndices[k] = corona.ringDisks[k]->type();
        diskPresence[diskIndices[k]] = true;
    }

    // count pairwise (adjacent) transitions (cyclic)
    for (size_t k = 0; k < diskIndices.size(); ++k) {
        size_t a = diskIndices[k];
        size_t b = diskIndices[(k + 1) % diskIndices.size()];
        ++getAdjacency(a, b);
    }
}

size_t& CoronaSignature::getAdjacency(size_t i, size_t j) {
    if (i > j) std::swap(i, j);
    size_t idx = (j * (j + 1)) / 2 + i;
    return adjacencyCounts[idx];
}

size_t CoronaSignature::getAdjacency(size_t i, size_t j) const {
    if (i > j) std::swap(i, j);
    size_t idx = (j * (j + 1)) / 2 + i;
    return adjacencyCounts[idx];
}

size_t CoronaSignature::baseType() const {
    return baseDiskType;
}

bool CoronaSignature::operator<(const CoronaSignature& other) const {
    if (baseDiskType != other.baseDiskType) return baseDiskType < other.baseDiskType;
    if (diskIndices.size() != other.diskIndices.size()) return diskIndices.size() < other.diskIndices.size();
    return adjacencyCounts < other.adjacencyCounts;
}

bool CoronaSignature::operator==(const CoronaSignature& other) const {
    return baseDiskType == other.baseDiskType &&
           diskIndices == other.diskIndices &&
           adjacencyCounts == other.adjacencyCounts;
}

// Reconstruct a corona from signature using operator table
bool CoronaSignature::testRadii(OperatorLookupTable& table) const {
    if (diskIndices.empty()) return false;

    // Build a base disk and first neighbor
    Disk base(0.0, 0.0, table.radii[baseDiskType], static_cast<size_t>(baseDiskType));
    std::vector<DiskPointer> packing;
    // place first disk to the right of base (deterministic)
    double x0 = (table.radii[baseDiskType].lower() + table.radii[baseDiskType].upper()) / 2.0
              + (table.radii[diskIndices[0]].lower() + table.radii[diskIndices[0]].upper()) / 2.0;
    packing.push_back(std::make_shared<Disk>( makeIntervalPoint(x0, 0.0), table.radii[diskIndices[0]], diskIndices[0] ));

    Corona testCorona(base, packing, table);

    // attempt to rebuild ring using the recorded sequence
    size_t frontIdx = 0;
    size_t backIdx = 0;

    for (size_t step = 1; step < diskIndices.size(); ++step) {
        bool useFront = testCorona.ringDisks.front()->precision() < testCorona.ringDisks.back()->precision();
        size_t nextType = diskIndices[ useFront ? frontIdx : backIdx ];

        auto placement = testCorona.calculatePlacement(nextType, std::nullopt);
        if (!placement.success) return false;

        // check intersections with existing packing
        for (const auto &dptr : packing) {
            if (placement.disk.intersects(*dptr)) return false;
        }

        packing.push_back(std::make_shared<Disk>(placement.disk));
        testCorona.addDisk(packing.back(), nextType);

        if (useFront) {
            frontIdx = (frontIdx == 0 ? diskIndices.size() - 1 : frontIdx - 1);
        } else {
            backIdx = (backIdx == diskIndices.size() - 1 ? 0 : backIdx + 1);
        }
    }

    return testCorona.isComplete();
}

// -------------------------------
// ConnectivityGraph
// -------------------------------

ConnectivityGraph::ConnectivityGraph(OperatorLookupTable& table)
    : precisionThreshold(table.radii.empty() ? BaseType(1e-6) : table.radii[0].lower() / BaseType(3.0)),
      overflowFlag(false),
      signaturesByBase(table.radii.size()),
      versionHistory(table.radii.size()),
      transitionCounts(table.radii.size(),
          std::vector<size_t>((table.radii.size() * table.radii.size() + table.radii.size()) / 2, 0)),
      adjacencyMatrix(table.radii.size(), std::vector<bool>(table.radii.size(), false))
{
    std::vector<DiskPointer> packing;
    // enumerate possible start disks and attempt to build coronas for each base type
    for (size_t baseType = 0; baseType < table.radii.size(); ++baseType) {
        Disk baseDisk(0.0, 0.0, table.radii[baseType], baseType);
        std::set<CoronaSignaturePointer> uniqueSignatures;

        for (size_t startType = 0; startType < table.radii.size(); ++startType) {
            // initial neighbor placed to the right
            double dx = (table.radii[baseType].lower() + table.radii[baseType].upper()) / 2.0
                      + (table.radii[startType].lower() + table.radii[startType].upper()) / 2.0;
            packing.push_back(std::make_shared<Disk>( makeIntervalPoint(dx, 0.0), table.radii[startType], startType ));

            Corona corona(baseDisk, packing, table);
            fillCorona(corona, packing, startType, uniqueSignatures);
            packing.clear();
        }
    }

    updateAdjacency();
}

void ConnectivityGraph::fillCorona(Corona& corona,
                                   std::vector<DiskPointer>& packing,
                                   size_t startIndex,
                                   std::set<CoronaSignaturePointer>& uniqueSignatures)
{
    if (corona.isComplete()) {
        auto signature = std::make_shared<CoronaSignature>(corona);
        if (uniqueSignatures.insert(signature).second) {
            signaturesByBase[corona.baseDisk().type()].push_back(signature);

            // update transitions
            const auto &indices = signature->diskIndices;
            for (size_t i = 0; i + 1 < indices.size(); ++i) {
                ++getTransition(corona.baseDisk().type(), indices[i], indices[i+1]);
            }
            if (!indices.empty()) {
                ++getTransition(corona.baseDisk().type(), indices.back(), indices.front());
            }
        }
        return;
    }

    // grow corona exploring all disk types >= startIndex (heuristic consistent with original)
    for (size_t diskType = startIndex; diskType < opTable.radii.size(); ++diskType) {
        auto placement = corona.calculatePlacement(diskType);
        if (!placement.success) continue;

        bool intersects = false;
        for (const auto &dptr : packing) {
            if (placement.disk.intersects(*dptr)) { intersects = true; break; }
        }
        if (intersects) continue;

        packing.push_back(std::make_shared<Disk>(placement.disk));
        corona.addDisk(packing.back(), diskType);

        fillCorona(corona, packing, startIndex, uniqueSignatures);

        corona.removeDisk();
        packing.pop_back();

        if (overflowFlag) return;
    }
}

size_t& ConnectivityGraph::getTransition(size_t base, size_t i, size_t j) {
    if (i > j) std::swap(i, j);
    size_t idx = (j * (j + 1)) / 2 + i;
    return transitionCounts[base][idx];
}

size_t ConnectivityGraph::getTransition(size_t base, size_t i, size_t j) const {
    if (i > j) std::swap(i, j);
    size_t idx = (j * (j + 1)) / 2 + i;
    return transitionCounts[base][idx];
}

size_t ConnectivityGraph::size() const {
    size_t s = 0;
    for (const auto &lst : signaturesByBase) s += lst.size();
    return s;
}

void ConnectivityGraph::display() const {
    std::cout << "ConnectivityGraph: bases = " << signaturesByBase.size()
              << " total signatures = " << size() << std::endl;
}

bool ConnectivityGraph::hasTriangle(size_t i, size_t j, size_t k) const {
    // naive check: if adjacencyMatrix indicates edges for these indices
    if (i >= adjacencyMatrix.size() || j >= adjacencyMatrix.size() || k >= adjacencyMatrix.size()) return false;
    return adjacencyMatrix[i][j] && adjacencyMatrix[j][k] && adjacencyMatrix[k][i];
}

void ConnectivityGraph::updateAdjacency() {
    // clear adjacency
    for (auto &row : adjacencyMatrix) std::fill(row.begin(), row.end(), false);

    // map triangular indices back to pairs and mark adjacency if transition > 0
    for (size_t b = 0; b < transitionCounts.size(); ++b) {
        const auto &vec = transitionCounts[b];
        for (size_t idx = 0; idx < vec.size(); ++idx) {
            if (vec[idx] == 0) continue;
            size_t i, j;
            triangularIndexToPair(idx, i, j);
            if (i < adjacencyMatrix.size() && j < adjacencyMatrix.size()) {
                adjacencyMatrix[i][j] = adjacencyMatrix[j][i] = true;
            }
        }
    }
}

bool ConnectivityGraph::isViable() const {
    // quick connectivity check on adjacency matrix (requires adjacency to be non-empty)
    if (adjacencyMatrix.empty()) return false;
    std::vector<bool> visited(adjacencyMatrix.size(), false);
    std::queue<size_t> q;
    q.push(0);
    visited[0] = true;
    while (!q.empty()) {
        size_t cur = q.front(); q.pop();
        for (size_t nb = 0; nb < adjacencyMatrix.size(); ++nb) {
            if (adjacencyMatrix[cur][nb] && !visited[nb]) {
                visited[nb] = true;
                q.push(nb);
            }
        }
    }
    return !std::any_of(visited.begin(), visited.end(), [](bool v){ return !v; });
}

void ConnectivityGraph::processInvalidTriangles(std::shared_ptr<std::vector<std::shared_ptr<SignatureList>>> /*diff*/) {
    // placeholder: real algorithm would pop invalid triangles and attempt to repair graph
    while (!invalidTriangles.empty()) invalidTriangles.pop();
}

bool ConnectivityGraph::restore() {
    // placeholder restore semantics
    // attempt to undo last version in versionHistory
    for (auto &vs : versionHistory) {
        if (!vs.empty()) {
            vs.pop();
        }
    }
    return true;
}

// -------------------------------
// Corona
// -------------------------------

Corona::Corona(const Disk& base, const std::vector<DiskPointer>& packing,
               OperatorLookupTable& table)
    : centralDisk(base),
      opTable(table)
{
    frontOperators.reserve(DEFAULT_OPERATORS);
    backOperators.reserve(DEFAULT_OPERATORS);

    // initialize ringDisks from packing
    for (const auto &dptr : packing) ringDisks.push_back(dptr);
    sortDisks();

    if (!ringDisks.empty()) {
        frontLeaf = { Interval(ringDisks.front()->x().lower() - centralDisk.x().lower(),
                               ringDisks.front()->x().upper() - centralDisk.x().upper()),
                      Interval(ringDisks.front()->y().lower() - centralDisk.y().lower(),
                               ringDisks.front()->y().upper() - centralDisk.y().upper()) };

        backLeaf  = { Interval(ringDisks.back()->x().lower() - centralDisk.x().lower(),
                               ringDisks.back()->x().upper() - centralDisk.x().upper()),
                      Interval(ringDisks.back()->y().lower() - centralDisk.y().lower(),
                               ringDisks.back()->y().upper() - centralDisk.y().upper()) };
    } else {
        frontLeaf = backLeaf = { Interval(0.0), Interval(0.0) };
    }
}

bool Corona::isContinuous() const {
    // continuous means each neighbor chain is touching next (non-gapped)
    if (ringDisks.size() < 2) return false;
    for (size_t i = 0; i + 1 < ringDisks.size(); ++i) {
        if (!ringDisks[i]->isTangent(*ringDisks[i + 1])) return false;
    }
    return true;
}

bool Corona::isComplete() const {
    if (ringDisks.size() < 3) return false;
    // ensure last touches first and every disk touches neighbors
    if (!ringDisks.front()->isTangent(*ringDisks.back())) return false;
    for (size_t i = 0; i < ringDisks.size(); ++i) {
        size_t next = (i + 1) % ringDisks.size();
        if (!ringDisks[i]->isTangent(*ringDisks[next])) return false;
    }
    return true;
}

Corona::Placement Corona::calculatePlacement(size_t diskType,
                                             const std::optional<ConnectivityGraph>& graph)
{
    // default: failure
    Placement out(false, Disk());

    if (ringDisks.empty()) return out;

    // Quick triangle check using connectivity graph if available
    if (graph) {
        if (!graph->hasTriangle(centralDisk.type(), ringDisks.front()->type(), diskType)) {
            return out;
        }
    }

    // decide building side (front/back) by precision heuristic
    bool useFront = ringDisks.front()->precision() < ringDisks.back()->precision();
    auto &ops = useFront ? frontOperators : backOperators;
    auto &leaf = useFront ? frontLeaf : backLeaf;

    // fetch operator from table (may throw or return identity if missing)
    SSORef ref = opTable.getOperator(centralDisk.type(),
                                     useFront ? ringDisks.front()->type() : ringDisks.back()->type(),
                                     diskType);
    ops.push_back(ref);

    // compose operator chain
    SpiralSimilarityOperator op = composeOperators(0, ops.size(), ops);
    if (useFront) {
        // convention: flip Y to place relative to front side
        op.flipY();
    }

    // apply operator to leaf vector (IntervalPair -> pair<Interval,Interval>)
    IntervalPair res = op * leaf; // assume operator * IntervalPair yields IntervalPair-like (first,second)
    // convert to absolute coordinates by adding central disk coords
    Interval cx = res.first + Interval(centralDisk.x().lower(), centralDisk.x().upper());
    Interval cy = res.second + Interval(centralDisk.y().lower(), centralDisk.y().upper());

    Disk placed;
    placed.center.x = cx;
    placed.center.y = cy;
    placed.radius = opTable.radii[diskType];
    placed.setType(diskType);

    // pop operator (we only trial-placed)
    ops.pop_back();

    out.success = true;
    out.disk = placed;
    return out;
}

SpiralSimilarityOperator Corona::composeOperators(size_t start, size_t end,
                                                  const std::vector<SSORef>& ops) const
{
    // Balanced binary composition to keep depth low
    size_t count = end - start;
    if (count == 0) return SpiralSimilarityOperator();
    if (count == 1) return ops[start];
    if (count == 2) return ops[start] * ops[start + 1];
    size_t mid = start + count / 2;
    return composeOperators(start, mid, ops) * composeOperators(mid, end, ops);
}

void Corona::addDisk(const DiskPointer& disk, size_t /*diskType*/) {
    ringDisks.push_back(disk);
    operationHistory.push(true);
    // keep ring sorted for deterministic order
    sortDisks();
}

void Corona::removeDisk() {
    if (ringDisks.empty()) return;
    ringDisks.pop_back();
    if (!operationHistory.empty()) operationHistory.pop();
}

void Corona::sortDisks() {
    if (ringDisks.size() <= 1) return;
    std::sort(ringDisks.begin(), ringDisks.end(),
        [this](const DiskPointer &a, const DiskPointer &b) {
            // choose representative coordinate (midpoint) for ordering
            double ax = (a->center.x.lower() + a->center.x.upper()) / 2.0 - (centralDisk.center.x.lower() + centralDisk.center.x.upper()) / 2.0;
            double ay = (a->center.y.lower() + a->center.y.upper()) / 2.0 - (centralDisk.center.y.lower() + centralDisk.center.y.upper()) / 2.0;
            double bx = (b->center.x.lower() + b->center.x.upper()) / 2.0 - (centralDisk.center.x.lower() + centralDisk.center.x.upper()) / 2.0;
            double by = (b->center.y.lower() + b->center.y.upper()) / 2.0 - (centralDisk.center.y.lower() + centralDisk.center.y.upper()) / 2.0;
            return std::atan2(ay, ax) < std::atan2(by, bx);
        });
}

const Disk& Corona::baseDisk() const {
    return centralDisk;
}

void Corona::display() const {
    std::cout << "Corona around base type " << centralDisk.type()
              << "  disks: " << ringDisks.size() << std::endl;
}

} // namespace diskpack