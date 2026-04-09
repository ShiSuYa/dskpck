#pragma once

#include <diskpack/geometry.h>
#include <diskpack/corona.h>

#include <set>
#include <random>
#include <optional>

namespace diskpack {

using RadiusType = BaseType;

struct LessNormCompare {
    bool operator()(const DiskPtr& a, const DiskPtr& b) const {
        return a->getNorm() < b->getNorm();
    }
};

/// ------------------------------------------------------------
/// Status of the packing generation process
/// ------------------------------------------------------------

enum class PackingStatus {

    complete,          /// packing successfully generated
    invalid,           /// no packing exists for these radii
    precision_error,   /// interval precision bound exceeded
    corona_error       /// corona assumption violated
};



/// ------------------------------------------------------------
/// Priority queue type
/// ------------------------------------------------------------
/// Returns disks closest to origin first

using DiskQueue =
    std::multiset<DiskPtr, LessNormCompare>;



std::ostream& operator<<(std::ostream& out, PackingStatus status);





/// ------------------------------------------------------------
/// Basic packing generator
/// ------------------------------------------------------------

class BasicGenerator {

protected:

    std::mt19937 rng;



    /// Radii of disk types
    std::vector<Interval> radii;



    /// Upper bound on interval width
    const RadiusType precision_upper_bound;



    /// Target radius of generated packing
    RadiusType packing_radius;



    /// Upper bound on disk count
    size_t size_upper_bound;



    /// Current packing state
    std::list<DiskPtr> packing;



    /// Queue of disks awaiting corona generation
    DiskQueue disk_queue;



    /// Lookup table for spiral operators
    SpiralOpCache op_cache;



    /// Disk frequency tracker
    std::vector<size_t> frequency_table;



    size_t max_ignored_radii;



    /// Radius covered by generated packing
    RadiusType generated_radius;



    /// --------------------------------------------------------
    /// Core recursive functions
    /// --------------------------------------------------------

    PackingStatus GapFill(Corona& corona);

    PackingStatus AdvancePacking();



    /// --------------------------------------------------------
    /// Utility functions
    /// --------------------------------------------------------

    void ShuffleIndexes(std::vector<size_t>& shuffle);

    void Push(Disk&& new_disk, size_t index);

    void Pop(size_t index);

    void SetGeneratedRadius(const Disk& furthest_disk);



    bool HasIntersection(const Disk& new_disk) const;

    bool IsInBounds(const Disk& disk) const;

    bool PackingSatisfiesConstraints() const;

    bool PackingIsLargeEnough() const;



public:

    std::optional<ConnectivityGraph> graph;



    BasicGenerator(
        const std::vector<Interval>& radii,
        const RadiusType& packing_radius,
        const RadiusType& precision_upper_bound,
        const size_t& size_upper_bound,
        const size_t& max_ignored_radii = 0
    );



    /// Start generation from scratch
    PackingStatus Generate(const size_t& central_disk_type);



    /// Reset generator state
    void Reset();



    /// Continue generation after parameter change
    PackingStatus Resume();



    /// Parameter setters
    void SetPackingRadius(const RadiusType& new_packing_radius);

    void SetSizeUpperBound(const size_t& new_size);

    void SetRadii(const std::vector<Interval>& radii_);



    /// Getters

    const RadiusType& GetGeneratedRadius();

    const RadiusType& GetRadius();

    const std::list<DiskPtr>& GetPacking();

};

} // namespace diskpack