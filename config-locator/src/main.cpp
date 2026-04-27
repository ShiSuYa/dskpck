#include <diskpack/search.h>
#include <diskpack/codec.h>

#include <boost/program_options.hpp>

#include <iomanip>
#include <fstream>
#include <iostream>
#include <chrono>
#include <thread>

using namespace diskpack;
namespace po = boost::program_options;

const size_t DEFAULT_SIZE_UPPER_BOUND = 25;
const BaseType DEFAULT_PACKING_RADIUS = 4;
const BaseType DEFAULT_PRECISION_UPPER_BOUND = 0.2;
const BaseType DEFAULT_LOWER_BOUND = 0.0000001;
const BaseType DEFAULT_UPPER_BOUND = 0.008;

int main(int argc, char *argv[]) {

    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    std::string output_file = "";

    size_t size_upper_bound, k;

    BaseType precision_upper_bound, packing_radius, lower_bound, upper_bound;

    RadiusRegion region{
      std::vector<Interval> {

        {0.5, 0.6},
        {0.6, 0.7},
        {0.7, 0.8},

        one, 
      }
    };

    try {
        po::options_description desc("Allowed options");

        desc.add_options()
            ("help,h", "Display usage information\n")
            
            ("lower-bound,l", 
                po::value<BaseType>()->default_value(DEFAULT_LOWER_BOUND), 
                "Minimum region width at which viability is checked")

            ("upper-bound,u", 
                po::value<BaseType>()->default_value(DEFAULT_UPPER_BOUND), 
                "Maximum region width to consider during the search")

            ("concurrency,k", 
                po::value<size_t>(), 
                "Number of worker threads to use")

            ("input,i", 
                po::value<std::string>(), 
                "Input JSON file describing the initial region")

            ("output,o", 
                po::value<std::string>(), 
                "Output JSON file for storing results");
    
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
                  std::cout << std::setprecision(7) << desc << "\n";

                  std::cout << "Example:\n";
                  std::cout << argv[0] 
                            << " -i input.json -n 100 -r 10\n\n";
    
                  return 0;
            }

            po::notify(vm);

            lower_bound = vm["lower-bound"].as<BaseType>();
            upper_bound = vm["upper-bound"].as<BaseType>();

            if (vm.count("concurrency")) {
                k = vm["concurrency"].as<size_t>();
            } else {
                k = std::thread::hardware_concurrency();
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
                  throw std::runtime_error(
                    "A single region must be specified"
                  );
                }

                region = regions[0].getIntervals();
              }

            if (vm.count("output")) {
                output_file = vm["output"].as<std::string>();
            }

      } catch (const po::error& e) {

        std::cerr << "Error: " << e.what() << "\n";
        return 1;

      } catch (const std::exception& e) {

        std::cerr << "Unhandled exception: " << e.what() << "\n";
        return 1;
    }

    std::vector<RadiusRegion> results;

    Searcher searcher{results, lower_bound, upper_bound};

    std::cerr << "Processing region:\t" 
              << EncodeRegionsJSON(std::vector<RadiusRegion> {region});

    std::cerr << "Threads:\t" 
              << std::thread::hardware_concurrency() 
              << "\n";

    auto t1 = high_resolution_clock::now();

    searcher.StartProcessing(region.getIntervals(), k);

    auto t2 = high_resolution_clock::now();
  
    auto ms_int = duration_cast<milliseconds>(t2 - t1);

    std::cerr << "Execution time:\t" 
              << ms_int.count()/60'000 
              << "m " 
              << (ms_int.count()/1000)%60 
              << "s\n";

    std::cerr << "Results size:\t" << results.size() << "\n";

    auto encoded = EncodeRegionsJSON(results);

    if (output_file == "") {

      std::cout << encoded << "\n";

    } else {

        std::ofstream out(output_file);

        if (!out.is_open()) {

            std::cerr << "Failed to open file: " + output_file + "\n";

            std::cout << encoded << "\n";

        } else {

            out << encoded << "\n";

        }

    }

    return 0;
}