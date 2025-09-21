#include <diskpack/search.h>
#include <diskpack/codec.h>

#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>

using namespace diskpack;
namespace po = boost::program_options;

namespace {
    constexpr size_t DEFAULT_SIZE_UPPER_BOUND      = 25;
    constexpr BaseType DEFAULT_PACKING_RADIUS      = 4;
    constexpr BaseType DEFAULT_PRECISION_UPPER_BOUND = 0.2;
    constexpr BaseType DEFAULT_LOWER_BOUND         = 0.0000001;
    constexpr BaseType DEFAULT_UPPER_BOUND         = 0.008;
}

int main(int argc, char* argv[]) {
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    std::string output_file;
    size_t k = std::thread::hardware_concurrency();
    BaseType lower_bound = DEFAULT_LOWER_BOUND;
    BaseType upper_bound = DEFAULT_UPPER_BOUND;

    RadiiRegion region{{
        {0.50, 0.65},
        {0.50, 0.65},
        {0.50, 0.65},
        one
    }};

    try {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "Show this help message")
            ("lower-bound,l", po::value<BaseType>()->default_value(DEFAULT_LOWER_BOUND),
                "Region width lower bound for viability check in the search")
            ("upper-bound,u", po::value<BaseType>()->default_value(DEFAULT_UPPER_BOUND),
                "Region width upper bound for viability check in the search")
            ("concurrency,k", po::value<size_t>(), "Number of threads in the thread pool")
            ("input,i", po::value<std::string>(), "Path to JSON file with the region search in")
            ("output,o", po::value<std::string>(), "Path to JSON file to store the output in");

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

        if (vm.count("help")) {
            std::cout << "Usage: " << argv[0] << "\n"
                      << std::setprecision(7) << desc << "\n"
                      << "Examples:\n"
                      << argv[0] << " -p 0.2 -n 200 -r 15 -l 0.001 -u 0.1\n\n"
                      << "The finder identifies all small subregions within a given region "
                         "which contain a radii set allowing compact disk packing.\n"
                      << "Use -n/-r to limit number of disks and covered region size.\n"
                      << "Use -p to set upper bound on interval precision.\n"
                      << "Regions narrower than lower_bound are added to results; "
                         "regions wider than upper_bound are skipped.\n";
            return 0;
        }

        po::notify(vm);

        lower_bound = vm["lower-bound"].as<BaseType>();
        upper_bound = vm["upper-bound"].as<BaseType>();

        if (vm.count("concurrency"))
            k = vm["concurrency"].as<size_t>();

        if (vm.count("input")) {
            std::ifstream file(vm["input"].as<std::string>());
            if (!file)
                throw std::runtime_error("Failed to open input file");

            std::vector<RadiiRegion> regions;
            DecodeRegionsJSON(file, regions);
            if (regions.size() != 1)
                throw std::runtime_error("Exactly one region should be provided");

            region = regions[0].GetIntervals();
        }

        if (vm.count("output"))
            output_file = vm["output"].as<std::string>();

    } catch (const po::error& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Unhandled exception: " << e.what() << "\n";
        return 1;
    }

    std::vector<RadiiRegion> results;
    Searcher searcher{results, lower_bound, upper_bound};

    std::cerr << "finder called on:\t" << EncodeRegionsJSON({region}) << "\n"
              << "hardware concurrency:\t" << std::thread::hardware_concurrency() << "\n";

    auto t1 = high_resolution_clock::now();
    searcher.StartProcessing(region.GetIntervals(), k);
    auto t2 = high_resolution_clock::now();

    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cerr << "duration:\t" << ms_int.count() / 60'000 << "m "
              << (ms_int.count() / 1000) % 60 << "s\n"
              << "results size:\t" << results.size() << "\n";

    const auto encoded = EncodeRegionsJSON(results);
    if (output_file.empty()) {
        std::cout << encoded << "\n";
    } else {
        std::ofstream out(output_file);
        if (!out)
            std::cerr << "Failed to open output file: " << output_file << "\n", std::cout << encoded << "\n";
        else
            out << encoded << "\n";
    }

    return 0;
}