#pragma once

#include <diskpack/geometry.h>
#include <list>

namespace diskpack {

    void WritePackingSVG(
      const std::string &filename, 
      const std::list<DiskPtr> &packing, 
      BaseType packing_radius
    );

    std::string EncodeRegionsJSON(
      const std::vector<RadiusRegion>& regions
    );

    void DecodeRegionsJSON(
      std::istream& data, 
      std::vector<RadiusRegion> &regions
    );
    
}