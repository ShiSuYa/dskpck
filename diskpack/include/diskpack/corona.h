#pragma once

#include <diskpack/geometry.h>

#include <queue>
#include <stack>
#include <set>
#include <list>
#include <deque>
#include <optional>

using DiskPointer = diskpack::DiskPtr;

namespace diskpack {

class Corona;
class CoronaSignature;
class ConnectivityGraph;

using CoronaSignaturePtr = std::shared_ptr<CoronaSignature>;
using CoronaSignatureList = std::list<CoronaSignaturePtr>;

bool CoronaSignaturePtrCompare(const CoronaSignaturePtr &a,
                               const CoronaSignaturePtr &b);

bool operator<(const CoronaSignaturePtr &a,
               const CoronaSignaturePtr &b);


/// ------------------------------------------------------------
/// ConnectivityGraph
///
/// Проверяет графовые инварианты, необходимые для существования
/// компактной упаковки.
/// ------------------------------------------------------------
class ConnectivityGraph {

public:

    const size_t MAX_SIGNATURES = 5000;

    BaseType PRECISION_THRESHOLD;

    using DiffStack = std::stack<std::shared_ptr<CoronaSignatureList>>;

private:

    bool broken_state = false;

    std::queue<std::tuple<size_t,size_t,size_t>> redundant_triangles;

    std::vector<CoronaSignatureList> signatures;

    std::vector<DiffStack> diffs;

    std::vector<std::vector<size_t>> transitions;

    std::vector<std::vector<bool>> adjacency_matrix;

private:

    void AddSignature(const CoronaSignature& signature);

    void RemoveSignature(const CoronaSignature& signature);

    void RemoveRedundantTriangles(
        std::shared_ptr<
            std::vector<std::shared_ptr<CoronaSignatureList>>
        > diff =
        std::shared_ptr<
            std::vector<std::shared_ptr<CoronaSignatureList>>
        >()
    );

    void UpdateEdges();

    size_t GetTransitionConst(size_t base,
                              size_t i,
                              size_t j) const;

    size_t& GetTransition(size_t base,
                          size_t i,
                          size_t j);

    void FillCorona(
        Corona& corona,
        std::list<DiskPointer> &packing,
        size_t start_index,
        std::set<
            CoronaSignaturePtr,
            decltype(&CoronaSignaturePtrCompare)
        > &unique_signatures
    );

public:

    size_t Size() const;

    void DisplaySignatures() const;

    ConnectivityGraph(SpiralOpCache &cache);

    bool HasOverflow() const;

    void Refine(SpiralOpCache &cache);

    bool Restore();

    bool IsViable() const;

    bool HasTriangle(size_t i,
                     size_t j,
                     size_t k) const;
};


/// ------------------------------------------------------------
/// Corona
///
/// Представляет корону дисков вокруг центрального.
/// ------------------------------------------------------------
class Corona {

    const size_t DEFAULT_OPERATOR_CAPACITY = 12;

    friend class CoronaSignature;
    friend class ConnectivityGraph;

private:

    const Disk &base_disk;

    std::deque<DiskPointer> corona_disks;

    std::vector<SpiralOpRef> operators_front;

    std::vector<SpiralOpRef> operators_back;

    IntervalPair leaf_front;

    IntervalPair leaf_back;

    SpiralOpCache &operator_cache;

    std::stack<bool> push_history;

private:

    SpiralOp computeOperatorsProduct(
        const size_t &begin,
        const size_t &end,
        const std::vector<SpiralOpRef> &ops
    ) const;

    void buildSortedCorona(
        const std::list<DiskPointer> &packing
    );

public:

    Corona(
        const Disk &base,
        const std::list<DiskPointer> &packing,
        SpiralOpCache &cache
    );

    bool isCompleted();

    bool isContinuous() const;

    bool peekNewDisk(
        Disk &new_disk,
        size_t index,
        const std::optional<ConnectivityGraph> &graph = std::nullopt
    );

    void PushDisk(
        const DiskPointer &disk,
        size_t index
    );

    void PopDisk();

    const Disk &getBase();

    void displaySignature();
};


/// ------------------------------------------------------------
/// CoronaSignature
///
/// Сигнатура короны: определяется количеством соседних пар
/// типов дисков.
/// ------------------------------------------------------------
class CoronaSignature {

private:

    const size_t radii_count;

    const size_t base_type;

    std::vector<size_t> transitions;

    size_t& GetTransition(size_t i,
                          size_t j);

    size_t GetTransitionConst(size_t i,
                              size_t j) const;

    friend class ConnectivityGraph;
    friend class Corona;

public:

    std::vector<bool> disk_presence;

    std::vector<size_t> specimen_indices;

    CoronaSignature(Corona &specimen);

    size_t getBase() const;

    bool TestRadii(SpiralOpCache& cache) const;

    bool operator<(const CoronaSignature &other) const;

    bool operator==(const CoronaSignature &other) const;
};

} // namespace diskpack