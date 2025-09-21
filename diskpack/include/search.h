#pragma once

#include <diskpack/generator.h>
#include <map>
#include <thread>
#include <vector>
#include <optional>

namespace diskpack {

class RegionAggregator {
    std::vector<size_t> componentSize_;
    std::vector<size_t> parent_;
    std::map<std::vector<Interval>, size_t, RadiiCompare> edgesMap_;
    std::vector<std::vector<Interval>> values_;

    size_t findParent(size_t x);
    void uniteSets(size_t x, size_t y);

public:
    RegionAggregator();
    void operator()(std::vector<RadiiRegion> &elements);
};

class RegionExplorer {
    std::vector<RadiiRegion>& resultsRef_;

    BaseType lowerBound_;   
    BaseType upperBound_;   

    bool expensiveCheck(const RadiiRegion& region);
    void processRegion(const RadiiRegion& region,
                       std::vector<RadiiRegion>& outRegions,
                       std::optional<ConnectivityGraph>& graph); 

public:
    RegionExplorer(std::vector<RadiiRegion> &results, BaseType lowerBound, BaseType upperBound);

    void startProcessing(const std::vector<Interval> &region, size_t concurrency = std::thread::hardware_concurrency());
};


}
