#include "geometry.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace diskpack {

bool Interval::Contains(BaseType value) const {
    return (low - EPSILON) <= value && value <= (high + EPSILON);
}

bool Interval::Overlaps(const Interval& other) const {
    return !(high < other.low - EPSILON || other.high < low - EPSILON);
}

BaseType Interval::Width() const {
    return high - low;
}

Interval Interval::Shifted(BaseType delta) const {
    return {low + delta, high + delta};
}

bool Interval::IsValid() const {
    return low <= high + EPSILON;
}

BaseType Disk::DistanceTo(const Disk& other) const {
    BaseType dx = center.x - other.center.x;
    BaseType dy = center.y - other.center.y;
    return std::sqrt(dx*dx + dy*dy);
}

BaseType Disk::SqDistanceTo(const Disk& other) const {
    BaseType dx = center.x - other.center.x;
    BaseType dy = center.y - other.center.y;
    return dx*dx + dy*dy;
}

bool Disk::Touches(const Disk& other, BaseType eps) const {
    BaseType dist = DistanceTo(other);
    BaseType sum = radius.high + other.radius.high;
    return std::abs(dist - sum) <= eps + EPSILON;
}

bool Disk::Intersects(const Disk& other) const {
    BaseType dist = DistanceTo(other);
    return dist < (radius.high + other.radius.high) - EPSILON;
}

bool Disk::Contains(const Point& p) const {
    BaseType dx = p.x - center.x;
    BaseType dy = p.y - center.y;
    return (dx*dx + dy*dy) <= (radius.high*radius.high) + EPSILON;
}

BaseType Disk::GetPrecision() const {
    return std::min({
        radius.Width(),
        std::abs(center.x - radius.low),
        std::abs(center.y - radius.low)
    });
}

Point DiskTransform::Apply(const Point& p) const {
    return {
        p.x * x.high - p.y * y.high,
        p.y * x.high + p.x * y.high
    };
}

DiskTransform DiskTransform::Compose(const DiskTransform& other) const {
    return {
        other.x.high * x.high - other.y.high * y.high,
        other.y.high * x.high + other.x.high * y.high
    };
}

DiskTransform DiskTransform::FromDisks(const Disk& base, const Disk& prev, const Disk& next) {
    const auto& br = base.radius;
    const auto& pr = prev.radius;
    const auto& nr = next.radius;

    Interval x;
    if (prev.type == base.type) {
        x = {0.5, 0.5};
    } 
    else if (next.type == base.type) {
        BaseType denom_high = 1.0 + pr.high/br.high;
        BaseType denom_low = 1.0 + pr.low/br.low;
        x = {2.0/denom_high, 2.0/denom_low};
    }
    else if (next.type == prev.type) {
        BaseType term_high = 2.0/(br.high/nr.high + 1.0);
        BaseType term_low = 2.0/(br.low/nr.low + 1.0);
        x = {1.0 - term_high, 1.0 - term_low};
    }
    else {
        BaseType ratio_high = pr.high/br.high;
        BaseType ratio_low = pr.low/br.low;
        
        BaseType first_term_high = 1.0/(1.0 + ratio_high);
        BaseType second_term_high = nr.high*(1.0 - ratio_high)/(br.high*pow(1.0 + ratio_high, 2));
        
        BaseType first_term_low = 1.0/(1.0 + ratio_low);
        BaseType second_term_low = nr.low*(1.0 - ratio_low)/(br.low*pow(1.0 + ratio_low, 2));
        
        x = {first_term_high + second_term_high, first_term_low + second_term_low};
    }

    Interval y;
    BaseType t_high = nr.high/(br.high + pr.high);
    BaseType t_low = nr.low/(br.low + pr.low);

    if (base.type == prev.type && next.type == prev.type) {
        y = {sqrt(3.0)/2.0, sqrt(3.0)/2.0};
    }
    else if (base.type == prev.type) {
        y = {
            sqrt(t_low * (1.0 + t_low)),
            sqrt(t_high * (1.0 + t_high))
        };
    }
    else {
        BaseType denom_high = 2.0 + pr.high/br.high + br.high/pr.high;
        BaseType denom_low = 2.0 + pr.low/br.low + br.low/pr.low;
        
        y = {
            2.0 * sqrt(t_low * (1.0 + t_low)) / denom_low,
            2.0 * sqrt(t_high * (1.0 + t_high)) / denom_high
        };
    }

    return {x, y};
}

DiskRegion::DiskRegion(const std::vector<Interval>& intervals) 
    : intervals(intervals) {
    if (intervals.empty()) {
        throw std::invalid_argument("Intervals list cannot be empty");
    }
}

const std::vector<Interval>& DiskRegion::GetIntervals() const {
    return intervals;
}

void DiskRegion::Split(std::vector<DiskRegion>& result, size_t parts, std::optional<size_t> idx) const {
    if (parts == 0) {
        throw std::invalid_argument("Number of parts must be positive");
    }

    result.clear();
    result.reserve(parts);

    auto intervals_copy = intervals;
    auto split_pos = idx.has_value() ? intervals_copy.begin() + *idx :
        std::max_element(intervals_copy.begin(), intervals_copy.end(),
            [](const Interval& a, const Interval& b) {
                return a.Width() < b.Width();
            });

    Interval original = *split_pos;
    BaseType step = original.Width() / parts;

    for (size_t i = 0; i < parts; ++i) {
        BaseType lower = original.low + i * step;
        BaseType upper = (i < parts - 1) ? lower + step : original.high;
        *split_pos = {lower, upper};
        result.emplace_back(intervals_copy);
    }
}

void DiskRegion::GridSplit(std::vector<DiskRegion>& result, size_t parts, size_t startIdx) const {
    if (startIdx == 0) {
        result.clear();
        size_t total = static_cast<size_t>(std::pow(parts, intervals.size() - 1));
        result.reserve(total);
    }

    if (startIdx >= intervals.size() - 1) {
        result.emplace_back(intervals);
        return;
    }

    std::vector<DiskRegion> temp_regions;
    Split(temp_regions, parts, startIdx);

    for (auto& region : temp_regions) {
        region.GridSplit(result, parts, startIdx + 1);
    }
}

bool DiskRegion::IsNarrow(BaseType threshold) const {
    return std::all_of(intervals.begin(), intervals.end(),
        [threshold](const Interval& i) {
            return i.Width() <= threshold + EPSILON;
        });
}

bool DiskRegion::IsWide(BaseType threshold) const {
    return std::any_of(intervals.begin(), intervals.end(),
        [threshold](const Interval& i) {
            return i.Width() > threshold - EPSILON;
        });
}

Interval DiskRegion::GetMinInterval() const {
    return *std::min_element(intervals.begin(), intervals.end(),
        [](const Interval& a, const Interval& b) {
            return a.low < b.low;
        });
}

Interval DiskRegion::GetMaxInterval() const {
    return *std::max_element(intervals.begin(), intervals.end(),
        [](const Interval& a, const Interval& b) {
            return a.high < b.high;
        });
}

DiskOperatorTable::DiskOperatorTable(const std::vector<Interval>& radii)
    : radii(radii),
      transforms(radii.size() * radii.size() * radii.size()),
      computed(radii.size() * radii.size() * radii.size(), false) {
    if (radii.empty()) {
        throw std::invalid_argument("Radii list cannot be empty");
    }
}

size_t DiskOperatorTable::GetIndex(size_t base, size_t prev, size_t next) const {
    return base + prev * radii.size() + next * radii.size() * radii.size();
}

const DiskTransform& DiskOperatorTable::GetTransform(size_t base, size_t prev, size_t next) {
    size_t idx = GetIndex(base, prev, next);
    if (!computed[idx]) {
        Disk b{{0, 0}, radii[base], base};
        Disk p{{0, 0}, radii[prev], prev};
        Disk n{{0, 0}, radii[next], next};
        transforms[idx] = DiskTransform::FromDisks(b, p, n);
        computed[idx] = true;
    }
    return transforms[idx];
}

namespace utils {

bool CompareByNorm(const Disk& a, const Disk& b) {
    BaseType na = a.center.x * a.center.x + a.center.y * a.center.y;
    BaseType nb = b.center.x * b.center.x + b.center.y * b.center.y;
    return (na < nb) || (std::abs(na - nb) < EPSILON && a.type < b.type);
}

bool CompareClockwise(const Disk& a, const Disk& b, const Point& center) {
    BaseType ax = a.center.x - center.x;
    BaseType ay = a.center.y - center.y;
    BaseType bx = b.center.x - center.x;
    BaseType by = b.center.y - center.y;

    if (ax > -EPSILON && std::abs(ay) <= EPSILON) return true;
    if (bx > -EPSILON && std::abs(by) <= EPSILON) return false;

    if (ay * by < -EPSILON) return ay > EPSILON;

    BaseType cross = bx * ay - by * ax;
    return cross > EPSILON || 
          (std::abs(cross) <= EPSILON && CompareByNorm(a, b));
}

}


}
