#pragma once

#include <diskpack/generator.h>

#include <map>
#include <thread>
#include <vector>
#include <optional>

namespace diskpack {

struct RadiusRegionCompare {
    bool operator()(const RadiusRegion& a, const RadiusRegion& b) const {
        return a.getIntervals() < b.getIntervals();
    }
};    


/// ------------------------------------------------------------
/// DSUFilter
/// ------------------------------------------------------------
///
/// Aggregates small regions into connected components.
/// Two regions are considered connected if they share
/// an adjacent side.
///
class DSUFilter {

private:

    std::vector<size_t> component_size;

    std::vector<size_t> parent;

    std::map<
        RadiusRegion,
        size_t,
        RadiusRegionCompare
    > edges;

    std::vector<std::vector<Interval>> vals;



    size_t Get(size_t x);

    void Unite(size_t x, size_t y);



public:

    DSUFilter();

    void operator()(std::vector<RadiusRegion>& elements);
};





/// ------------------------------------------------------------
/// Searcher
/// ------------------------------------------------------------
///
/// Performs exhaustive search in the radius parameter space.
/// Splits regions recursively and runs packing generation
/// checks for each subregion.
///
class Searcher {

private:

    std::vector<RadiusRegion>& results;



    /// If region diameter ≤ lower_bound → accept
    RadiusType lower_bound;



    /// If region diameter ≥ upper_bound → skip expensive check
    RadiusType upper_bound;



    bool ExpensiveCheck(const RadiusRegion& region);



    void ProcessRegion(
        const RadiusRegion& region,
        std::vector<RadiusRegion>& r,
        std::optional<ConnectivityGraph>& graph
    );



public:

    Searcher(
        std::vector<RadiusRegion>& results,
        RadiusType lower_bound,
        RadiusType upper_bound
    );



    /// Start parallel search
    void StartProcessing(
        std::vector<Interval> region,
        size_t k = std::thread::hardware_concurrency()
    );
};

} // namespace diskpack