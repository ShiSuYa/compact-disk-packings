#ifndef DISKPACK_CODEC_H
#define DISKPACK_CODEC_H

#include <string>
#include <vector>
#include <list>
#include <iosfwd>
#include "geometry.h"

namespace diskpack {

// SVG экспорт упаковки дисков
void exportToSVG(const std::string& filename, 
                const std::list<Disk>& disks,
                double boundaryRadius);

// Кодирование/декодирование областей
struct RegionData {
    std::vector<Point> centers;
    std::vector<Interval> radii;
};

std::string encodeRegions(const std::vector<DiskRegion>& regions);
bool decodeRegions(const std::string& json, RegionData& output);
bool decodeRegions(std::istream& stream, RegionData& output);

} // namespace diskpack

#endif // DISKPACK_CODEC_H