#include "codec.h"
#include <fstream>
#include <nlohmann/json.hpp>

namespace diskpack {

namespace {

const std::vector<std::string> DISK_COLORS = {
    "#ff8080a0", "#aa0087a0", "#ff0066a0", "#803dc1a0", "#000000a0"
};

void writeSvgHeader(std::ofstream& svg, double radius) {
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
        << "<svg viewBox=\"" << -radius << " " << -radius << " " 
        << 2 * radius << " " << 2 * radius 
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n"
        << "  <rect x=\"" << -radius << "\" y=\"" << -radius << "\" "
        << "width=\"" << 2 * radius << "\" height=\"" << 2 * radius << "\" "
        << "fill=\"white\"/>\n"
        << "  <g transform=\"scale(1, -1)\">\n";
}

}

void exportToSVG(const std::string& filename, 
                const std::list<Disk>& disks,
                double boundaryRadius) {
    
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    writeSvgHeader(file, boundaryRadius);

    for (const auto& disk : disks) {
        file << "  <circle"
             << " cx=\"" << disk.center.x << "\""
             << " cy=\"" << disk.center.y << "\""
             << " r=\"" << disk.radius.high << "\""
             << " fill=\"" << DISK_COLORS[disk.type % DISK_COLORS.size()] << "\"/>\n";
    }

    file << "  </g>\n</svg>\n";
}

std::string encodeRegions(const std::vector<DiskRegion>& regions) {
    nlohmann::json result;
    
    for (const auto& region : regions) {
        nlohmann::json intervals;
        for (const auto& ival : region.GetIntervals()) {
            intervals.push_back({ival.low, ival.high});
        }
        result.push_back(intervals);
    }

    return result.dump();
}

bool decodeRegions(const std::string& json, RegionData& output) {
    try {
        auto data = nlohmann::json::parse(json);
        return decodeRegionsInternal(data, output);
    } catch (const nlohmann::json::exception& e) {
        return false;
    }
}

bool decodeRegions(std::istream& stream, RegionData& output) {
    try {
        auto data = nlohmann::json::parse(stream);
        return decodeRegionsInternal(data, output);
    } catch (const nlohmann::json::exception& e) {
        return false;
    }
}

namespace {

bool decodeRegionsInternal(const nlohmann::json& data, RegionData& output) {
    if (!data.is_array()) return false;

    output.centers.clear();
    output.radii.clear();

    for (const auto& item : data) {
        if (!item.is_array() || item.size() != 2) continue;
        
        try {
            output.radii.emplace_back(
                item[0].get<double>(),
                item[1].get<double>()
            );
        } catch (...) {
            return false;
        }
    }

    return !output.radii.empty();
}

}


}
