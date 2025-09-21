#pragma once

#include <diskpack/corona.h>
#include <set>
#include <random>

namespace diskpack {

enum class GenerationResult {
    Success,          
    Impossible,       
    PrecisionLimit,   
    CoronaViolation   
};

using DiskPriorityQueue = std::multiset<DiskPointer, decltype(&LessNormCompare)>;

std::ostream& operator<<(std::ostream& out, GenerationResult result);

class PackingBuilderBase {
protected:
    std::mt19937 rngEngine;

    std::vector<Interval> diskRadii;          
    const BaseType maxPrecisionWidth;         
    BaseType targetPackingRadius;             
    size_t maxDiskCount;                      
    std::list<DiskPointer> currentPacking;    
    DiskPriorityQueue pendingDisks;           
    OperatorLookupTable operatorTable;        
    std::vector<size_t> diskUsageCounter;     
    size_t ignoreRadiusLimit;                 
    BaseType achievedPackingRadius;           

    GenerationResult fillCorona(Corona& corona);

    GenerationResult stepForward();

    void shuffleIndices(std::vector<size_t>& order);

    void pushDisk(Disk&& disk, size_t typeIndex);

    void popDisk(size_t typeIndex);

    void updateAchievedRadius(const Disk& furthest);

    bool intersectsExisting(const Disk& candidate) const;

    bool fitsInsideBounds(const Disk& disk) const;

    bool checkConstraints() const;

    bool largeEnough() const;

public:
    std::optional<ConnectivityGraph> connectivityMap;

    PackingBuilderBase(const std::vector<Interval>& radii,
                       const BaseType& targetRadius,
                       const BaseType& precisionLimit,
                       const size_t& diskCountLimit,
                       const size_t& maxIgnored = 0);

    GenerationResult generate(const size_t& centerDiskType);

    void reset();

    GenerationResult resume();

    void setTargetRadius(const BaseType& newRadius);

    void setDiskCountLimit(const size_t& newLimit);

    void setRadii(const std::vector<Interval>& newRadii);

    const BaseType& getAchievedRadius();

    const BaseType& getTargetRadius();

    const std::list<DiskPointer>& getPacking() const;
};


}
