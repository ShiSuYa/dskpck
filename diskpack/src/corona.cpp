#include <diskpack/corona.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>

namespace diskpack {

/// ============================================================
/// Corona
/// ============================================================

SpiralOp
Corona::computeOperatorsProduct(
    const size_t &begin,
    const size_t &end,
    const std::vector<SpiralOpRef> &operators
) const {

    switch (end - begin) {

    case 0:
        return SpiralOp();

    case 1:
        return operators[begin];

    case 2:
        return operators[begin].get() *
               operators[begin + 1].get();

    case 3:
        return operators[begin].get() *
               (operators[begin + 1].get() *
                operators[begin + 2].get());

    case 4:
        return (operators[begin].get() *
                operators[begin + 1].get()) *
               (operators[begin + 2].get() *
                operators[begin + 3].get());
    }

    size_t mid = (begin + end) / 2;

    return computeOperatorsProduct(begin, mid, operators) *
           computeOperatorsProduct(mid, end, operators);
}

bool Corona::isCompleted() {

    assert(!corona_disks.empty());

    Interval cross_product =
        (corona_disks.back()->getCenterX() - base_disk.getCenterX()) *
        (corona_disks.front()->getCenterY() - base_disk.getCenterY()) -
        (corona_disks.front()->getCenterX() - base_disk.getCenterX()) *
        (corona_disks.back()->getCenterY() - base_disk.getCenterY());

    if (!(corona_disks.size() > 2) ||
        !corona_disks.back()->tangent(*corona_disks.front()) ||
        !cergt(cross_product, 0.0L)) {
        return false;
    }

    bool use_front =
        corona_disks.front()->precision() <
        corona_disks.back()->precision();

    Disk new_disk;

    const Disk &old_disk =
        *(!use_front ? corona_disks.front() : corona_disks.back());

    peekNewDisk(new_disk, old_disk.getType());

    bool check_intersect =
        !empty(intersect(new_disk.getCenterX(), old_disk.getCenterX())) &&
        !empty(intersect(new_disk.getCenterY(), old_disk.getCenterY()));

    return check_intersect;
}

void Corona::buildSortedCorona(const std::list<DiskPointer> &packing) {

    assert(corona_disks.empty());

    for (auto &disk : packing) {
        if (base_disk.tangent(*disk)) {
            corona_disks.push_back(disk);
        }
    }

    std::sort(
        corona_disks.begin(),
        corona_disks.end(),
        ClockwiseDiskCompare{base_disk}
    );

    for (size_t i = 0; i + 1 < corona_disks.size(); ++i) {
        if (!corona_disks[i]->tangent(*corona_disks[i + 1])) {
            std::rotate(
                corona_disks.begin(),
                corona_disks.begin() + i + 1,
                corona_disks.end()
            );
            return;
        }
    }
}

bool Corona::isContinuous() const {

    for (size_t i = 0; i + 1 < corona_disks.size(); ++i) {
        if (!corona_disks[i]->tangent(*corona_disks[i + 1])) {
            return false;
        }
    }
    return true;
}

Corona::Corona(
    const Disk &base,
    const std::list<DiskPointer> &packing,
    SpiralOpCache &cache
)
    : base_disk(base),
      operator_cache(cache) {

    operators_back.reserve(DEFAULT_OPERATOR_CAPACITY);
    operators_front.reserve(DEFAULT_OPERATOR_CAPACITY);

    buildSortedCorona(packing);

    assert(!corona_disks.empty());

    leaf_front =
        IntervalPair{
            corona_disks.front()->getCenterX() - base_disk.getCenterX(),
            corona_disks.front()->getCenterY() - base_disk.getCenterY()
        };

    leaf_back =
        IntervalPair{
            corona_disks.back()->getCenterX() - base_disk.getCenterX(),
            corona_disks.back()->getCenterY() - base_disk.getCenterY()
        };
}

bool Corona::peekNewDisk(
    Disk &new_disk,
    size_t index,
    const std::optional<ConnectivityGraph> &graph
) {

    bool use_front =
        corona_disks.front()->precision() <
        corona_disks.back()->precision();

    if (graph.has_value()) {
        if (!graph->HasTriangle(
                base_disk.getType(),
                use_front ?
                corona_disks.front()->getType() :
                corona_disks.back()->getType(),
                index)) {
            return false;
        }
    }

    auto &operators =
        use_front ? operators_front : operators_back;

    auto &start_leaf =
        use_front ? leaf_front : leaf_back;

    operators.push_back(
        operator_cache(
            base_disk.getType(),
            use_front ?
            corona_disks.front()->getType() :
            corona_disks.back()->getType(),
            index
        )
    );

    SpiralOp op =
        computeOperatorsProduct(0, operators.size(), operators);

    op.y *= (use_front ? -1 : 1);

    IntervalPair new_center = op * start_leaf;

    new_disk = Disk(
        new_center.first + base_disk.getCenterX(),
        new_center.second + base_disk.getCenterY(),
        operator_cache.radii[index],
        index
    );

    operators.pop_back();

    return true;
}

void Corona::PushDisk(const DiskPointer &disk, size_t index) {

    bool use_front =
        corona_disks.front()->precision() <
        corona_disks.back()->precision();

    push_history.push(use_front);

    auto &operators =
        use_front ? operators_front : operators_back;

    operators.push_back(
        operator_cache(
            base_disk.getType(),
            use_front ?
            corona_disks.front()->getType() :
            corona_disks.back()->getType(),
            index
        )
    );

    if (use_front)
        corona_disks.push_front(disk);
    else
        corona_disks.push_back(disk);
}

void Corona::PopDisk() {

    bool use_front = push_history.top();
    push_history.pop();

    if (use_front) {
        operators_front.pop_back();
        corona_disks.pop_front();
    } else {
        operators_back.pop_back();
        corona_disks.pop_back();
    }
}

const Disk &Corona::getBase() {
    return base_disk;
}

void Corona::displaySignature() {

    std::cerr << "signature:\n";

    CoronaSignature signature(*this);

    std::cerr << "base: " << signature.getBase() << "\n";

    for (size_t i = 0; i < operator_cache.radii.size(); ++i) {
        for (size_t j = 0; j < operator_cache.radii.size(); ++j) {
            std::cerr << signature.GetTransitionConst(i, j) << " ";
        }
        std::cerr << "\n";
    }

    for (auto &i : signature.specimen_indices) {
        std::cerr << i << " ";
    }

    std::cerr << "\n";
}

/// ============================================================
/// CoronaSignature (DTSP)
/// ============================================================

size_t& CoronaSignature::GetTransition(size_t i, size_t j) {
    return transitions[i > j ? (i * i + i) / 2 + j
                             : (j * j + j) / 2 + i];
}

size_t CoronaSignature::GetTransitionConst(size_t i, size_t j) const {
    return transitions[i > j ? (i * i + i) / 2 + j
                             : (j * j + j) / 2 + i];
}

bool CoronaSignature::operator<(const CoronaSignature &other) const {
    return transitions < other.transitions;
}

bool CoronaSignature::operator==(const CoronaSignature &other) const {
    return transitions == other.transitions;
}

CoronaSignature::CoronaSignature(Corona &specimen)
    : radii_count(specimen.operator_cache.radii.size()),
      base_type(specimen.base_disk.getType()),
      transitions((radii_count * radii_count + radii_count) / 2),
      disk_presence(radii_count, false),
      specimen_indices(specimen.corona_disks.size()) {

    if (!specimen.isCompleted()) {
        throw std::runtime_error("CoronaSignature requires completed corona");
    }

    for (size_t i = 0; i < specimen_indices.size(); ++i) {
        specimen_indices[i] =
            specimen.corona_disks[i]->getType();
        disk_presence[specimen_indices[i]] = true;
    }

    for (size_t i = 0; i + 1 < specimen_indices.size(); ++i) {
        ++GetTransition(specimen_indices[i],
                        specimen_indices[i + 1]);
    }

    ++GetTransition(specimen_indices.back(),
                    specimen_indices.front());
}

size_t CoronaSignature::getBase() const {
    return base_type;
}

bool CoronaSignaturePtrCompare(const CoronaSignaturePtr &a,
                               const CoronaSignaturePtr &b) {
    return *a < *b;
}

bool CoronaSignature::TestRadii(SpiralOpCache &cache) const {

    auto &radii = cache.radii;

    Disk base_disk(0, 0, radii[base_type], base_type);

    std::list<DiskPointer> packing;

    packing.push_back(
        std::make_shared<Disk>(
            radii[base_type] + radii[specimen_indices[0]],
            0,
            radii[specimen_indices[0]],
            specimen_indices[0]
        )
    );

    Corona test(base_disk, packing, cache);

    size_t cur_front = 0, cur_back = 0;

    for (size_t i = 1; i < specimen_indices.size(); ++i) {

        bool use_front =
            test.corona_disks.front()->precision() <
            test.corona_disks.back()->precision();

        if (use_front)
            cur_front =
                (cur_front == 0 ?
                 specimen_indices.size() - 1 :
                 cur_front - 1);
        else
            cur_back =
                (cur_back == specimen_indices.size() - 1 ?
                 0 :
                 cur_back + 1);

        size_t next =
            specimen_indices[use_front ? cur_front : cur_back];

        Disk new_disk;

        test.peekNewDisk(new_disk, next);

        if (std::any_of(
                packing.begin(),
                packing.end(),
                [&new_disk](const DiskPointer &d) {
                    return new_disk.intersects(*d);
                }))
        {
            return false;
        }

        packing.push_back(
            std::make_shared<Disk>(new_disk)
        );

        test.PushDisk(packing.back(), next);
    }

    return test.isCompleted();
}

/// ------------------------------------------------------------
/// ConnectivityGraph (FULL DTSP - COMPLETE)
/// ------------------------------------------------------------

void ConnectivityGraph::FillCorona(
    Corona& corona,
    std::list<DiskPointer>& packing,
    size_t start_index,
    std::set<CoronaSignaturePtr, decltype(&CoronaSignaturePtrCompare)>& unique
) {
    auto& radii = corona.operator_cache.radii;

    if (corona.isCompleted()) {
        auto signature = std::make_shared<CoronaSignature>(corona);

        if (unique.insert(signature).second) {
            signatures[corona.getBase().getType()].push_back(signature);
            AddSignature(*signature);
        }
        return;
    }

    Disk new_disk;

    for (size_t i = start_index; i < radii.size(); ++i) {

        corona.peekNewDisk(new_disk, i);

        if (std::any_of(packing.begin(), packing.end(),
            [&new_disk](const DiskPointer& d) {
                return new_disk.intersects(*d);
            }))
        {
            continue;
        }

        if (new_disk.precision() > PRECISION_THRESHOLD) {
            broken_state = true;
            return;
        }

        packing.push_back(std::make_shared<Disk>(std::move(new_disk)));
        corona.PushDisk(packing.back(), i);

        FillCorona(corona, packing, start_index, unique);

        corona.PopDisk();
        packing.pop_back();

        if (HasOverflow()) return;
    }
}


/// ------------------------------------------------------------
/// Constructor
/// ------------------------------------------------------------
ConnectivityGraph::ConnectivityGraph(SpiralOpCache& cache)
    : diffs(cache.radii.size()),
      signatures(cache.radii.size()),
      transitions(cache.radii.size(),
          std::vector<size_t>(
              (cache.radii.size() * cache.radii.size() + cache.radii.size()) / 2,
              0)),
      adjacency_matrix(cache.radii.size(),
          std::vector<bool>(cache.radii.size(), false))
{
    std::list<DiskPointer> packing;
    auto& radii = cache.radii;

    PRECISION_THRESHOLD =
        std::min_element(radii.begin(), radii.end(),
            [](const Interval& a, const Interval& b) {
                return a.lower() < b.lower();
            })->lower() / 3.0L;

    for (size_t base = 0; base < radii.size(); ++base) {

        Disk base_disk(0, 0, radii[base], base);

        std::set<CoronaSignaturePtr, decltype(&CoronaSignaturePtrCompare)>
            unique(CoronaSignaturePtrCompare);

        for (size_t start = 0; start < radii.size(); ++start) {

            packing.push_back(
                std::make_shared<Disk>(
                    radii[base] + radii[start], 0,
                    radii[start], start
                )
            );

            Corona corona(base_disk, packing, cache);

            FillCorona(corona, packing, start, unique);

            packing.clear();
        }
    }

    for (size_t base = 0; base < adjacency_matrix.size(); ++base) {
        for (size_t i = 0; i < adjacency_matrix.size(); ++i) {
            for (size_t j = i; j < adjacency_matrix.size(); ++j) {
                if (GetTransitionConst(base, i, j) == 0) {
                    redundant_triangles.push({base, i, j});
                }
            }
        }
    }

    if (HasOverflow()) return;

    RemoveRedundantTriangles();
    UpdateEdges();
}


/// ------------------------------------------------------------
size_t ConnectivityGraph::Size() const {
    size_t res = 0;
    for (size_t i = 0; i < signatures.size(); ++i)
        res += signatures[i].size();
    return res;
}


/// ------------------------------------------------------------
bool ConnectivityGraph::HasOverflow() const {
    return broken_state ||
        std::any_of(signatures.begin(), signatures.end(),
            [this](const CoronaSignatureList& lst) {
                return lst.size() > MAX_SIGNATURES;
            });
}


/// ------------------------------------------------------------
bool ConnectivityGraph::HasTriangle(size_t i, size_t j, size_t k) const {
    return GetTransitionConst(i, j, k) &&
           GetTransitionConst(k, i, j) &&
           GetTransitionConst(j, k, i);
}


/// ------------------------------------------------------------
size_t& ConnectivityGraph::GetTransition(size_t base, size_t i, size_t j) {
    return transitions[base][
        i > j ? (i * i + i) / 2 + j
              : (j * j + j) / 2 + i
    ];
}

size_t ConnectivityGraph::GetTransitionConst(size_t base, size_t i, size_t j) const {
    return transitions[base][
        i > j ? (i * i + i) / 2 + j
              : (j * j + j) / 2 + i
    ];
}


/// ------------------------------------------------------------
void ConnectivityGraph::AddSignature(const CoronaSignature& sig) {
    for (size_t x = 0; x < sig.transitions.size(); ++x) {
        transitions[sig.getBase()][x] += sig.transitions[x];
    }
}


/// ------------------------------------------------------------
void ConnectivityGraph::RemoveSignature(const CoronaSignature& sig) {
    for (size_t i = 0; i < adjacency_matrix.size(); ++i) {
        for (size_t j = i; j < adjacency_matrix.size(); ++j) {

            GetTransition(sig.getBase(), i, j) -=
                sig.GetTransitionConst(i, j);

            if (sig.GetTransitionConst(i, j) > 0 &&
                GetTransitionConst(sig.getBase(), i, j) == 0)
            {
                redundant_triangles.push({i, j, sig.getBase()});
            }
        }
    }
}


/// ------------------------------------------------------------
void ConnectivityGraph::RemoveRedundantTriangles(
    std::shared_ptr<std::vector<std::shared_ptr<CoronaSignatureList>>> diff)
{
    std::vector<CoronaSignatureList> removed(signatures.size());

    while (!redundant_triangles.empty()) {

        auto [i,j,k] = redundant_triangles.front();
        redundant_triangles.pop();

        std::vector<size_t> triple{i,j,k};

        for (size_t r = 0; r < 3; ++r) {

            std::rotate(triple.begin(), triple.begin()+1, triple.end());

            if (GetTransitionConst(triple[0], triple[1], triple[2]) == 0)
                continue;

            for (auto it = signatures[triple[0]].begin();
                 it != signatures[triple[0]].end();)
            {
                auto sig = *it;

                if (sig->GetTransitionConst(triple[1], triple[2]) > 0) {

                    removed[triple[0]].push_back(sig);
                    it = signatures[triple[0]].erase(it);

                    RemoveSignature(*sig);

                } else {
                    ++it;
                }
            }
        }
    }

    if (diff) {
        for (size_t i = 0; i < signatures.size(); ++i) {
            diff->at(i)->splice(diff->at(i)->end(), std::move(removed[i]));
        }
    }
}


/// ------------------------------------------------------------
void ConnectivityGraph::UpdateEdges() {

    for (size_t i = 0; i < adjacency_matrix.size(); ++i) {
        for (size_t j = i; j < adjacency_matrix.size(); ++j) {

            adjacency_matrix[i][j] = false;

            for (size_t x = 0; x < adjacency_matrix.size(); ++x) {

                adjacency_matrix[i][j] =
                    adjacency_matrix[i][j] ||
                    (GetTransitionConst(i,j,x) > 0 &&
                     GetTransitionConst(j,i,x) > 0);

                adjacency_matrix[j][i] = adjacency_matrix[i][j];
            }
        }
    }
}


/// ------------------------------------------------------------
bool ConnectivityGraph::Restore() {

    if (std::any_of(diffs.begin(), diffs.end(),
        [](const DiffStack& s){ return s.empty(); }))
        return false;

    for (size_t base = 0; base < diffs.size(); ++base) {

        auto diff = diffs[base].top();
        diffs[base].pop();

        for (auto sig : *diff)
            AddSignature(*sig);

        signatures[base].splice(signatures[base].end(), std::move(*diff));
    }

    UpdateEdges();
    return true;
}


/// ------------------------------------------------------------
void ConnectivityGraph::Refine(SpiralOpCache& cache) {

    auto invalid =
        std::make_shared<
            std::vector<std::shared_ptr<CoronaSignatureList>>
        >(diffs.size());

    for (auto& l : *invalid)
        l = std::make_shared<CoronaSignatureList>();

    for (size_t base = 0; base < cache.radii.size(); ++base) {

        for (auto it = signatures[base].begin();
             it != signatures[base].end();)
        {
            auto sig = *it;

            if (sig->TestRadii(cache)) {
                ++it;
                continue;
            }

            invalid->at(base)->push_back(sig);
            it = signatures[base].erase(it);

            RemoveSignature(*sig);
        }
    }

    RemoveRedundantTriangles(invalid);

    for (size_t base = 0; base < diffs.size(); ++base) {
        diffs[base].push(invalid->at(base));
    }

    UpdateEdges();
}

/// ------------------------------------------------------------
bool ConnectivityGraph::IsViable() const {

    /// ------------------------------------------------------------
    /// Condition 1: Graph is connected (у тебя уже было — оставляем)
    /// ------------------------------------------------------------
    {
        std::vector<bool> vis(adjacency_matrix.size(), false);
        std::queue<size_t> q;

        q.push(0);
        vis[0] = true;

        while (!q.empty()) {
            auto v = q.front(); q.pop();

            for (size_t u = 0; u < adjacency_matrix.size(); ++u) {
                if (!vis[u] && adjacency_matrix[v][u]) {
                    vis[u] = true;
                    q.push(u);
                }
            }
        }

        if (std::any_of(vis.begin(), vis.end(),
                        [](bool b){ return !b; }))
            return false;
    }

    /// ------------------------------------------------------------
    /// Condition 2
    /// ------------------------------------------------------------
    {
        for (size_t removed = 0; removed < adjacency_matrix.size(); ++removed) {

            std::vector<size_t> components(adjacency_matrix.size(), 0);
            components[removed] = adjacency_matrix.size();

            size_t current_component = 0;

            for (size_t v = 0; v < adjacency_matrix.size(); ++v) {

                if (components[v]) continue;

                ++current_component;

                std::queue<size_t> q;
                q.push(v);
                components[v] = current_component;

                while (!q.empty()) {
                    size_t x = q.front(); q.pop();

                    for (size_t u = 0; u < adjacency_matrix.size(); ++u) {
                        if (components[u] || !adjacency_matrix[x][u]) continue;

                        components[u] = current_component;
                        q.push(u);
                    }
                }
            }

            if (!std::any_of(
                    signatures[removed].begin(),
                    signatures[removed].end(),
                    [current_component, this, &components]
                    (const CoronaSignaturePtr &sig)
                    {
                        for (size_t comp = 1; comp <= current_component; ++comp) {

                            bool connects = false;

                            for (size_t v = 0; v < adjacency_matrix.size(); ++v) {
                                connects = connects ||
                                    (sig->disk_presence[v] &&
                                     components[v] == comp);
                            }

                            if (!connects) return false;
                        }

                        return true;
                    }))
            {
                return false;
            }
        }
    }

    /// ------------------------------------------------------------
    /// Condition 3
    /// ------------------------------------------------------------
    {
        size_t N = adjacency_matrix.size();

        std::vector<size_t> parent(N * N * N);
        std::vector<size_t> comp_size(parent.size(), 1);

        std::iota(parent.begin(), parent.end(), 0);

        auto get_index = [&](size_t i, size_t j, size_t k) {
            std::vector<size_t> ijk{i, j, k};
            std::sort(ijk.begin(), ijk.end());
            return ijk[0] + N * ijk[1] + N * N * ijk[2];
        };

        std::function<size_t(size_t)> find = [&](size_t x) {
            return parent[x] == x ? x : parent[x] = find(parent[x]);
        };

        auto unite = [&](size_t a, size_t b) {
            a = find(a);
            b = find(b);
            if (a != b) {
                if (comp_size[a] < comp_size[b]) std::swap(a, b);
                parent[b] = a;
                comp_size[a] += comp_size[b];
            }
        };

        for (size_t base = 0; base < N; ++base) {

            for (auto &sig : signatures[base]) {

                size_t i0 = sig->specimen_indices[0];
                size_t j0 = sig->specimen_indices[1];

                size_t first = get_index(base, i0, j0);

                for (size_t i = 0; i < N; ++i) {
                    for (size_t j = i; j < N; ++j) {

                        if (sig->GetTransitionConst(i, j) > 0) {
                            unite(first, get_index(base, i, j));
                        }
                    }
                }
            }
        }

        std::vector<std::vector<bool>> comp_has(N * N * N,
            std::vector<bool>(N, false));

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i; j < N; ++j) {
                for (size_t k = j; k < N; ++k) {

                    size_t comp = find(get_index(i, j, k));

                    comp_has[comp][i] = true;
                    comp_has[comp][j] = true;
                    comp_has[comp][k] = true;
                }
            }
        }

        if (!std::any_of(
                comp_has.begin(),
                comp_has.end(),
                [](const std::vector<bool> &v) {
                    return std::all_of(v.begin(), v.end(),
                                       [](bool x){ return x; });
                }))
        {
            return false;
        }
    }

    return true;
}

} // namespace diskpack