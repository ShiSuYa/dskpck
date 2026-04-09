#pragma once

#include <diskpack/geometry.h>
#include <diskpack/generator.h>

#include <list>
#include <vector>
#include <string>
#include <istream>

namespace diskpack {

/// Export disk packing as SVG image
void WritePackingSVG(
    const std::string& filename,
    const std::list<DiskPtr>& packing,
    RadiusType packing_radius
);

/// Encode radius regions to JSON string
std::string EncodeRegionsJSON(
    const std::vector<RadiusRegion>& regions
);

/// Decode radius regions from JSON stream
void DecodeRegionsJSON(
    std::istream& data,
    std::vector<RadiusRegion>& regions
);

} // namespace diskpack