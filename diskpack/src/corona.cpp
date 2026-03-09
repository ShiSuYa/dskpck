#include <diskpack/corona.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>

namespace diskpack {

/// ------------------------------------------------------------
/// Helper: multiply a range of spiral operators
/// ------------------------------------------------------------

SpiralOp
Corona::computeOperatorsProduct(
    const size_t& begin,
    const size_t& end,
    const std::vector<SpiralOpRef>& operators
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






/// ------------------------------------------------------------
/// Check if corona is closed
/// ------------------------------------------------------------

bool Corona::isCompleted() {

    assert(!corona.empty());

    Interval cross_product =
        (corona.back()->getCenterX() - base.getCenterX()) *
        (corona.front()->getCenterY() - base.getCenterY()) -

        (corona.front()->getCenterX() - base.getCenterX()) *
        (corona.back()->getCenterY() - base.getCenterY());


    if (!(corona.size() > 2) ||
        !corona.back()->tangent(*corona.front()) ||
        !cergt(cross_product, 0.0L))
    {
        return false;
    }

    bool use_front =
        corona.front()->precision() <
        corona.back()->precision();

    Disk new_disk;

    const Disk& old_disk =
        *(!use_front ? corona.front() : corona.back());

    peekNewDisk(new_disk, old_disk.getType());

    bool check_intersect =
        !empty(intersect(new_disk.getCenterX(), old_disk.getCenterX())) &&
        !empty(intersect(new_disk.getCenterY(), old_disk.getCenterY()));

    return check_intersect;
}






/// ------------------------------------------------------------
/// Build ordered corona from packing
/// ------------------------------------------------------------

void Corona::buildSortedCorona(const std::list<DiskPtr>& packing) {

    assert(corona.empty());

    for (auto& disk : packing) {

        if (base.tangent(*disk)) {
            corona.push_back(disk);
        }
    }

    std::sort(
        corona.begin(),
        corona.end(),
        ClockwiseDiskCompare{base}
    );

    for (size_t i = 0; i + 1 < corona.size(); ++i) {

        if (!corona[i]->tangent(*corona[i + 1])) {

            std::rotate(
                corona.begin(),
                corona.begin() + i + 1,
                corona.end()
            );

            return;
        }
    }
}






/// ------------------------------------------------------------
/// Corona continuity check
/// ------------------------------------------------------------

bool Corona::isContinuous() const {

    for (size_t i = 0; i + 1 < corona.size(); ++i) {

        if (!corona[i]->tangent(*corona[i + 1])) {
            return false;
        }
    }

    return true;
}






/// ------------------------------------------------------------
/// Constructor
/// ------------------------------------------------------------

Corona::Corona(
    const Disk& base_disk,
    const std::list<DiskPtr>& packing,
    SpiralOpCache& cache
)
    :
    base(base_disk),
    op_cache(cache)
{

    operators_back.reserve(DEFAULT_OPERATORS_SIZE);
    operators_front.reserve(DEFAULT_OPERATORS_SIZE);

    buildSortedCorona(packing);

    assert(!corona.empty());

    leaf_front =
        IntervalPair{
            corona.front()->getCenterX() - base.getCenterX(),
            corona.front()->getCenterY() - base.getCenterY()
        };

    leaf_back =
        IntervalPair{
            corona.back()->getCenterX() - base.getCenterX(),
            corona.back()->getCenterY() - base.getCenterY()
        };
}






/// ------------------------------------------------------------
/// Predict new disk without inserting
/// ------------------------------------------------------------

bool Corona::peekNewDisk(
    Disk& new_disk,
    size_t index,
    const std::optional<ConnectivityGraph>& graph
) {

    bool use_front =
        corona.front()->precision() <
        corona.back()->precision();


    if (graph.has_value()) {

        if (!graph->HasTriangle(
                base.getType(),
                use_front ?
                corona.front()->getType() :
                corona.back()->getType(),
                index))
        {
            return false;
        }
    }

    auto& operators =
        use_front ? operators_front : operators_back;

    auto& start_leaf =
        use_front ? leaf_front : leaf_back;


    operators.push_back(
        op_cache(
            base.getType(),
            use_front ?
            corona.front()->getType() :
            corona.back()->getType(),
            index
        )
    );


    SpiralOp op =
        computeOperatorsProduct(
            0,
            operators.size(),
            operators
        );


    op.y *= (use_front ? -1 : 1);


    IntervalPair new_center =
        op * start_leaf;


    new_disk = Disk(
        new_center.first + base.getCenterX(),
        new_center.second + base.getCenterY(),
        op_cache.radii[index],
        index
    );

    operators.pop_back();

    return true;
}






/// ------------------------------------------------------------
/// Push disk to corona
/// ------------------------------------------------------------

void Corona::push(const DiskPtr& disk, size_t index) {

    bool use_front =
        corona.front()->precision() <
        corona.back()->precision();

    push_history.push(use_front);

    auto& operators =
        use_front ? operators_front : operators_back;


    operators.push_back(
        op_cache(
            base.getType(),
            use_front ?
            corona.front()->getType() :
            corona.back()->getType(),
            index
        )
    );


    if (use_front)
        corona.push_front(disk);
    else
        corona.push_back(disk);
}






/// ------------------------------------------------------------
/// Pop disk from corona
/// ------------------------------------------------------------

void Corona::pop() {

    bool use_front = push_history.top();
    push_history.pop();

    if (use_front) {

        operators_front.pop_back();
        corona.pop_front();

    } else {

        operators_back.pop_back();
        corona.pop_back();
    }
}






/// ------------------------------------------------------------
/// Access base disk
/// ------------------------------------------------------------

const Disk& Corona::getBase() {

    return base;
}






/// ------------------------------------------------------------
/// Debug print
/// ------------------------------------------------------------

void Corona::displaySignature() {

    std::cerr << "signature:\n";

    CoronaSignature signature(*this);

    std::cerr << "base: "
              << signature.base
              << "\n";

    for (size_t i = 0; i < op_cache.radii.size(); ++i) {

        for (size_t j = 0; j < op_cache.radii.size(); ++j) {

            std::cerr
                << signature.getTransitions(i, j)
                << " ";
        }

        std::cerr << "\n";
    }

    for (auto& i : signature.specimen_indexes) {
        std::cerr << i << " ";
    }

    std::cerr << "\n";
}

} // namespace diskpack
