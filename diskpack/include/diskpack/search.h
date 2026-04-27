#pragma once

#include <diskpack/generator.h>

#include <map>
#include <thread>


namespace diskpack {

    class DSUFilter {     // DSUFilter merges adjacent regions into connected components. Each output region is a union of connected regions.

        std::vector<size_t> component_size;

        std::vector<size_t> parent;

        std::map<std::vector<Interval>, size_t, RadiiCompare> edges;

        std::vector<std::vector<Interval>> vals;

        size_t Get(size_t x);

        void Unite(size_t x, size_t y);

    public:

        DSUFilter();

        void operator()(std::vector<RadiusRegion> &elements);
    };
    

    class Searcher {      // Searches for radii that allow compact packings by splitting the region. Recursively refines subregions and stops when they become small enough.

        std::vector<RadiusRegion>& results;

        BaseType lower_bound;     // Regions smaller than this are accepted

        BaseType upper_bound;     // Regions larger than this are skipped

        bool ExpensiveCheck(const RadiusRegion& region);

        void ProcessRegion(
          const RadiusRegion& region, 
          std::vector<RadiusRegion>& r, 
          std::optional<ConnectivityGraph> &graph
        );     // Checks whether a region may contain valid radii
    
    public:

        Searcher(
          std::vector<RadiusRegion> &results, 
          BaseType lower_bound, 
          BaseType upper_bound
        );

        void StartProcessing(
          std::vector<Interval> region, 
          size_t k = std::thread::hardware_concurrency()
        );     // Starts the parallel search

    };
    
}