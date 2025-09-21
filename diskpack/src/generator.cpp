#include <diskpack/generator.h>
#include <iostream>
#include <numeric>  // для std::iota
#include <algorithm> // для std::shuffle, std::any_of, std::count
#include <cmath> // для sqrt

namespace diskpack {

// Используемое приватное поле генератора случайных чисел
static std::random_device random_device_source;

PackingBuilderBase::PackingBuilderBase(const std::vector<Interval>& radii_,
                                       const BaseType& targetRadius_,
                                       const BaseType& precisionLimit_,
                                       const size_t& maxDiskCount_,
                                       const size_t& maxIgnoredRadii_)
    : diskRadii(radii_),
      targetPackingRadius(targetRadius_),
      maxPrecisionWidth(precisionLimit_),
      maxDiskCount(maxDiskCount_),
      ignoreRadiusLimit(maxIgnoredRadii_),
      pendingDisks(LessNormCompare),
      diskUsageCounter(radii_.size(), 0),
      operatorTable(radii_),
      achievedPackingRadius(0),
      rngEngine(random_device_source()),
      connectivityMap(operatorTable)
{
    if (connectivityMap->HasOverflow()) {
        connectivityMap.reset();
    }
}

GenerationResult PackingBuilderBase::generate(const size_t& centerDiskType) {
    if (diskRadii.empty()) {
        throw std::runtime_error("PackingBuilderBase::generate() called with empty radii vector");
    }
    if (connectivityMap.has_value() && !connectivityMap->IsViable()) {
        std::cerr << "Connectivity map is not viable!\n";
        return GenerationResult::Impossible;
    }

    reset();

    // Начинаем упаковку с центрального диска
    pushDisk(Disk(zero, zero, diskRadii[centerDiskType], centerDiskType), centerDiskType);

    // Добавляем второй диск разных радиусов по очереди
    for (size_t i = 0; i < diskRadii.size(); ++i) {
        if (i == centerDiskType && diskRadii.size() > 1) {
            continue;
        }

        pushDisk(Disk(diskRadii[centerDiskType] + diskRadii[i], zero, diskRadii[i], i), i);

        auto status = stepForward();
        if (status != GenerationResult::Impossible) {
            return status;
        }

        popDisk(i);
    }

    popDisk(centerDiskType);
    return GenerationResult::Impossible;
}

GenerationResult PackingBuilderBase::resume() {
    if (currentPacking.size() < 2) {
        return GenerationResult::Impossible;
    }
    return stepForward();
}

void PackingBuilderBase::reset() {
    pendingDisks.clear();
    currentPacking.clear();
    achievedPackingRadius = 0;
    std::fill(diskUsageCounter.begin(), diskUsageCounter.end(), 0);
}

void PackingBuilderBase::shuffleIndices(std::vector<size_t>& order) {
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rngEngine);
}

GenerationResult PackingBuilderBase::stepForward() {
    if (pendingDisks.empty()) {
        return GenerationResult::Impossible;
    }

    auto baseDisk = pendingDisks.extract(pendingDisks.begin()).value();

    if (baseDisk->precision() > maxPrecisionWidth) {
        updateAchievedRadius(*baseDisk);
        return GenerationResult::PrecisionLimit;
    }

    if (!fitsInsideBounds(*baseDisk) || largeEnough()) {
        pendingDisks.insert(baseDisk);

        if (!checkConstraints()) {
            return GenerationResult::Impossible;
        }

        updateAchievedRadius(*baseDisk);
        return GenerationResult::Success;
    }

    Corona corona(*baseDisk, currentPacking, operatorTable);
    if (!corona.IsContinuous()) {
        updateAchievedRadius(*baseDisk);
        return GenerationResult::CoronaViolation;
    }

    auto status = fillCorona(corona);
    if (status == GenerationResult::Impossible) {
        pendingDisks.insert(baseDisk);
    }

    return status;
}

GenerationResult PackingBuilderBase::fillCorona(Corona& corona) {
    if (corona.IsCompleted()) {
        return stepForward();
    }

    std::vector<size_t> order(diskRadii.size());
    shuffleIndices(order);

    Disk candidateDisk;

    for (size_t i = 0; i < diskRadii.size(); ++i) {
        if (!corona.PeekNewDisk(candidateDisk, order[i], connectivityMap)) {
            continue;
        }

        if (intersectsExisting(candidateDisk)) {
            continue;
        }

        pushDisk(std::move(candidateDisk), order[i]);
        corona.Push(currentPacking.back(), order[i]);

        auto status = fillCorona(corona);
        if (status != GenerationResult::Impossible) {
            return status;
        }

        corona.Pop();
        popDisk(order[i]);
    }
    return GenerationResult::Impossible;
}

void PackingBuilderBase::pushDisk(Disk&& disk, size_t typeIndex) {
    currentPacking.push_back(std::make_shared<Disk>(std::move(disk)));
    pendingDisks.insert(currentPacking.back());
    ++diskUsageCounter[typeIndex];
}

void PackingBuilderBase::popDisk(size_t typeIndex) {
    --diskUsageCounter[typeIndex];
    pendingDisks.erase(currentPacking.back());
    currentPacking.pop_back();
}

bool PackingBuilderBase::intersectsExisting(const Disk& candidate) const {
    return std::any_of(currentPacking.begin(), currentPacking.end(),
        [&candidate](const DiskPointer& dptr) {
            return candidate.intersects(*dptr);
        });
}

bool PackingBuilderBase::fitsInsideBounds(const Disk& disk) const {
    BaseType dist = disk.get_norm();
    BaseType radius_margin = targetPackingRadius - disk.get_radius();
    return cerle(dist, radius_margin * radius_margin);
}

bool PackingBuilderBase::checkConstraints() const {
    return std::count(diskUsageCounter.begin(), diskUsageCounter.end(), 0) <= ignoreRadiusLimit;
}

bool PackingBuilderBase::largeEnough() const {
    return currentPacking.size() >= maxDiskCount;
}

void PackingBuilderBase::updateAchievedRadius(const Disk& furthest) {
    BaseType dist = sqrt(furthest.get_norm()) + furthest.get_radius();
    achievedPackingRadius = median(dist);
}

void PackingBuilderBase::setTargetRadius(const BaseType& newRadius) {
    targetPackingRadius = newRadius;
}

void PackingBuilderBase::setDiskCountLimit(const size_t& newLimit) {
    maxDiskCount = newLimit;
}

void PackingBuilderBase::setRadii(const std::vector<Interval>& newRadii) {
    reset();
    diskRadii = newRadii;
    operatorTable = OperatorLookupTable(newRadii);
    diskUsageCounter.resize(newRadii.size(), 0);
}

const BaseType& PackingBuilderBase::getAchievedRadius() {
    return achievedPackingRadius;
}

const BaseType& PackingBuilderBase::getTargetRadius() {
    return targetPackingRadius;
}

const std::list<DiskPointer>& PackingBuilderBase::getPacking() const {
    return currentPacking;
}

std::ostream& operator<<(std::ostream& out, GenerationResult result) {
    switch (result) {
    case GenerationResult::Success:
        return out << "Success";
    case GenerationResult::Impossible:
        return out << "Impossible";
    case GenerationResult::PrecisionLimit:
        return out << "PrecisionLimit";
    case GenerationResult::CoronaViolation:
        return out << "CoronaViolation";
    }
    return out;
}

} // namespace diskpack