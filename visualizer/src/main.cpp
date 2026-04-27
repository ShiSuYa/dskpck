#include <diskpack/constants.h>
#include <diskpack/generator.h>
#include <diskpack/codec.h>

#include <boost/program_options.hpp>

#include <chrono>
#include <iostream>
#include <fstream>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

using namespace diskpack;
namespace po = boost::program_options;

const size_t DEFAULT_SIZE_UPPER_BOUND = 200;
const BaseType DEFAULT_PACKING_RADIUS = 5;
const BaseType DEFAULT_PRECISION_UPPER_BOUND = 0.5;
const std::string DEFAULT_OUTPUT_FILE = "../assets/default.svg";

std::vector<Interval> radii;

int main(int argc, char *argv[]) {

  bool use_custom_values = false;

  size_t size_upper_bound;
  size_t central_disk_type;

  BaseType packing_radius;
  BaseType precision_upper_bound;

  std::string output_file;

  try {

    po::options_description desc("Allowed options");

    desc.add_options()
        ("help,h", "Show this help message\n")

        ("region-size,r", 
            po::value<BaseType>()->default_value(DEFAULT_PACKING_RADIUS), 
            "Upper limit on the region size")

        ("number-of-disks,n", 
            po::value<size_t>()->default_value(DEFAULT_SIZE_UPPER_BOUND),
            "Upper limit on the number of disks")

        ("output,o", 
            po::value<std::string>()->default_value(DEFAULT_OUTPUT_FILE), 
            "SVG file to store the generated packing")
        
        ("precision,p", 
            po::value<BaseType>()->default_value(DEFAULT_PRECISION_UPPER_BOUND), 
            "Upper limit on coordinate precision")

        ("central-disk,c", 
            po::value<size_t>(), 
            "Central disk type")

        ("i2", 
            po::value<size_t>(), 
            "Use the i-th radius allowing packing with radii 1 and r")

        ("i3", 
            po::value<size_t>(), 
            "Use the i-th pair (r, s) allowing packing with radii 1, r, s")

        ("input,i", 
            po::value<std::string>(), 
            "JSON file with the region");

        po::positional_options_description p;

        po::variables_map vm;

        po::store(
          po::command_line_parser(argc, argv)
            .options(desc)
            .positional(p)
            .run(), 
          vm
        );

        if (vm.count("help")) {

              std::cout << "Usage of " << argv[0] << "\n";
              std::cout << desc << "\n";

              return 0;
        }

        po::notify(vm);

        size_upper_bound = vm["number-of-disks"].as<size_t>();
        precision_upper_bound = vm["precision"].as<BaseType>();
        packing_radius = vm["region-size"].as<BaseType>();
        output_file = vm["output"].as<std::string>();

        if ((vm.count("i2") && vm.count("i3")) || 
            (vm.count("input") && vm.count("i3")) || 
            (vm.count("i2") && vm.count("input"))  ||
            (!vm.count("i2") && !vm.count("i3") && !vm.count("input"))) {

            throw std::runtime_error(
              "Provide exactly one of the flags: i2, i3, or i"
            );
        }

        if (vm.count("i2")) {
          auto i2 = vm["i2"].as<size_t>() - 1;
          if (i2 < 0 || i2 >= 9) {
            throw std::runtime_error("i2 out of range");
          }
          radii = std::vector<Interval>{one, two__disk_radii[i2]};
        }

        if (vm.count("i3")) {
          auto i3 = vm["i3"].as<size_t>() - 1;
          if (i3 < 0 || i3 >= 164) {
            throw std::runtime_error("i3 out of range");
          }
          radii = std::vector<Interval>{one, three_disk_radii[i3].first, three_disk_radii[i3].second};
        }

        if (vm.count("input")) {

          std::vector<RadiusRegion> regions;

          std::string filename = vm["input"].as<std::string>();

          std::ifstream file(filename);

          if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
          }

          DecodeRegionsJSON(file, regions);

          if (regions.size() != 1) {
            throw std::runtime_error("A single region must be provided");
          }

          radii = regions[0].getIntervals();
        }

        if (vm.count("central-disk")) {

          central_disk_type = vm["central-disk"].as<size_t>();

        } else {

          central_disk_type = std::distance(
            radii.begin(), 
            std::min_element(
              radii.begin(), 
              radii.end(), 
              [](const Interval &a, const Interval &b) {
                return a.lower() < b.lower();
              }
            )
          );
        }
        
  } catch (const po::error& e) {

    std::cerr << "Error: " << e.what() << "\n";
    return 1;

  } catch (const std::exception& e) {

    std::cerr << "Unhandled exception: " << e.what() << "\n";
    return 1;
  }
  
  std::sort(
    radii.begin(), 
    radii.end(), 
    [](const Interval& a, const Interval& b) {
      return cerlt(a, b);
    }
  );

  std::cerr << "Processing region: " 
            << EncodeRegionsJSON(std::vector<RadiusRegion> {{radii}});

  BasicGenerator generator{
    radii, 
    packing_radius, 
    precision_upper_bound, 
    size_upper_bound
  };

  auto t1 = high_resolution_clock::now();

  auto status = generator.Generate(central_disk_type);

  auto t2 = high_resolution_clock::now();

  auto ms_int = duration_cast<milliseconds>(t2 - t1);
  
  std::cout << "Status:\t" << status << "\n";
  std::cout << "Execution time:\t" << ms_int.count() << "ms\n";

  if (status != PackingStatus::invalid) {

    WritePackingSVG(
      output_file, 
      generator.GetPacking(), 
      generator.GetGeneratedRadius() + 1);
  }

  return 0;
}