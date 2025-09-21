#include <diskpack/constants.h>
#include <diskpack/generator.h>
#include <diskpack/codec.h>

#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

using Clock      = std::chrono::high_resolution_clock;
using Millis     = std::chrono::milliseconds;
using diskpack::BaseType;
using diskpack::Interval;
using diskpack::RadiiRegion;
using diskpack::BasicGenerator;
using diskpack::PackingStatus;

namespace po = boost::program_options;

// Значения по умолчанию
constexpr size_t    SIZE_LIMIT_DEFAULT       = 200;
constexpr BaseType  RADIUS_LIMIT_DEFAULT     = 5;
constexpr BaseType  PRECISION_LIMIT_DEFAULT  = 0.5;
const std::string   DEFAULT_OUTPUT_PATH      = "../images/default.svg";

static std::vector<Interval> radiiSet;

int main(int argc, char* argv[]) {
    size_t    diskCountLimit   = SIZE_LIMIT_DEFAULT;
    size_t    centralDiskIndex = 0;
    BaseType  regionRadius     = RADIUS_LIMIT_DEFAULT;
    BaseType  precisionLimit   = PRECISION_LIMIT_DEFAULT;
    std::string outputPath     = DEFAULT_OUTPUT_PATH;

    try {
        // Определяем и описываем аргументы командной строки
        po::options_description opts("Параметры запуска");
        opts.add_options()
            ("help,h", "Показать справку")
            ("region-size,r", po::value<BaseType>()->default_value(RADIUS_LIMIT_DEFAULT),
                "Максимальный радиус области покрытия")
            ("number-of-disks,n", po::value<size_t>()->default_value(SIZE_LIMIT_DEFAULT),
                "Максимальное количество дисков")
            ("output,o", po::value<std::string>()->default_value(DEFAULT_OUTPUT_PATH),
                "Имя .svg файла для сохранения результата")
            ("precision,p", po::value<BaseType>()->default_value(PRECISION_LIMIT_DEFAULT),
                "Максимальная допустимая точность координат")
            ("central-disk,c", po::value<size_t>(), 
                "Тип центрального диска (индекс радиуса в наборе)")
            ("i2", po::value<size_t>(), "Использовать набор (1, r)")
            ("i3", po::value<size_t>(), "Использовать набор (1, r, s)")
            ("input,i", po::value<std::string>(), "JSON с описанием области поиска");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, opts), vm);

        if (vm.count("help")) {
            std::cout << "Использование: " << argv[0] << "\n" << opts << "\n";
            return 0;
        }
        po::notify(vm);

        // Читаем параметры
        diskCountLimit = vm["number-of-disks"].as<size_t>();
        precisionLimit = vm["precision"].as<BaseType>();
        regionRadius   = vm["region-size"].as<BaseType>();
        outputPath     = vm["output"].as<std::string>();

        // Проверка на корректность выбора источника радиусов
        const bool i2Set    = vm.count("i2");
        const bool i3Set    = vm.count("i3");
        const bool inputSet = vm.count("input");

        if ((i2Set + i3Set + inputSet) != 1) {
            throw std::runtime_error("Должен быть указан ровно один из {i2, i3, input}");
        }

        if (i2Set) {
            size_t idx = vm["i2"].as<size_t>() - 1;
            if (idx >= 9) throw std::runtime_error("i2: 0 < i <= 9");
            radiiSet = { diskpack::one, diskpack::two_radii[idx] };
        }
        if (i3Set) {
            size_t idx = vm["i3"].as<size_t>() - 1;
            if (idx >= 164) throw std::runtime_error("i3: 0 < i <= 164");
            radiiSet = { diskpack::one, diskpack::three_radii[idx].first, diskpack::three_radii[idx].second };
        }
        if (inputSet) {
            std::ifstream inFile(vm["input"].as<std::string>());
            if (!inFile) throw std::runtime_error("Не удалось открыть входной файл");
            std::vector<RadiiRegion> regions;
            diskpack::DecodeRegionsJSON(inFile, regions);
            if (regions.size() != 1) throw std::runtime_error("Должен быть ровно один регион");
            radiiSet = regions[0].GetIntervals();
        }

        if (vm.count("central-disk")) {
            centralDiskIndex = vm["central-disk"].as<size_t>();
        } else {
            centralDiskIndex = std::distance(
                radiiSet.begin(),
                std::min_element(radiiSet.begin(), radiiSet.end(),
                    [](const Interval& a, const Interval& b) {
                        return a.lower() < b.lower();
                    })
            );
        }

    } catch (const std::exception& ex) {
        std::cerr << "Ошибка: " << ex.what() << "\n";
        return 1;
    }

    // Сортируем радиусы
    std::sort(radiiSet.begin(), radiiSet.end(), diskpack::cerlt);
    std::cerr << "Visualizer data: " 
              << diskpack::EncodeRegionsJSON({ {radiiSet} });

    // Генерация упаковки
    BasicGenerator generator(radiiSet, regionRadius, precisionLimit, diskCountLimit);

    auto start = Clock::now();
    auto status = generator.Generate(centralDiskIndex);
    auto finish = Clock::now();

    auto elapsed = std::chrono::duration_cast<Millis>(finish - start);

    std::cout << "status:\t"   << status << "\n"
              << "time:\t"     << elapsed.count() << " ms\n";

    if (status != PackingStatus::invalid) {
        diskpack::WritePackingSVG(outputPath, generator.GetPacking(), generator.GetGeneratedRadius() + 1);
    }

    return 0;
}