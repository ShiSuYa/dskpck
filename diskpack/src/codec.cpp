#include <diskpack/codec.h>

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

const std::vector<std::string> type_to_color = {
    "#ff8080a0", 
    "#aa0087a0", 
    "#ff0066a0", 
    "#803dc1a0", 
    "#000000a0",
};

namespace diskpack {

void WritePackingSVG(
  const std::string &filename,
  const std::list<DiskPtr> &packing,
  BaseType packing_radius)  
{
  std::ofstream svgFile(filename);
  
  int width = 1000, height = 1000;

  if (!svgFile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  svgFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";

  svgFile << "<svg viewBox=\"" << -packing_radius << " " << -packing_radius
          << " " << 2 * packing_radius << " " << 2 * packing_radius
          << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";

  svgFile << "  <rect x=\"" << -packing_radius << "\" y=\"" << -packing_radius
          << "\" "
          << "width=\"" << 2 * packing_radius << "\" height=\""
          << 2 * packing_radius << "\" "
          << "fill=\"white\"/>\n";

  svgFile << "  <g transform=\"scale(1, -1)\">\n";

  for (const auto &disk : packing) {
    svgFile << "  <circle"
            << " cx=\"" << median(disk->getCenterX()) 
            << "\" cy=\"" << median(disk->getCenterY()) 
            << "\" r=\"" << median(disk->getRadius()) 
            << "\" fill=\""
            << type_to_color[disk->getType()] 
            << "\"/>\n";
  }

  svgFile << "  </g>\n";
  svgFile << "</svg>\n";

  svgFile.close();

  std::cerr << "SVG file created successfully: " 
            << filename << std::endl;
}

std::string EncodeRegionsJSON(
  const std::vector<RadiusRegion> &regions) 
{
  std::string result;

  for (const auto &region : regions) {
    json json_line = json::array();

    for (const auto &interval : region.getIntervals()) {
      json_line.push_back(
        json::array({
          interval.lower(), 
          interval.upper()
        })
      );
    }

    result += json_line.dump() + '\n';

  }

  return result;

}

void DecodeRegionsJSON(
  std::istream &data, 
  std::vector<RadiusRegion> &regions) 
{
  regions.clear();

  std::string region;

  while (std::getline(data, region)) {

    if (region.empty()) {
      continue;
    }

    json json_line;

    try {

      json_line = json::parse(region);

    } catch (const json::parse_error &e) {

      throw std::runtime_error(
        "JSON parse error: " + 
        std::string(e.what()));
    }

    if (!json_line.is_array()) {

      throw std::runtime_error(
        "Each JSON line must be an array"
      );
    }

    std::vector<Interval> current_line;

    for (const auto &json_point : json_line) {

      if (!json_point.is_array() || 
          json_point.size() != 2) 
      {
        throw std::runtime_error(
            "Invalid point format - expected array of two numbers");
      }

      current_line.emplace_back(
        json_point[0].get<long double>(),
        json_point[1].get<long double>()
      );
    }

    regions.emplace_back(current_line);
  }
}

}