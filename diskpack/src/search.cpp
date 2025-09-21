#include <diskpack/codec.h>
#include <diskpack/search.h>

#include <atomic>
#include <iostream>
#include <string>
#include <thread>
#include <algorithm>

namespace diskpack {

/// RegionAggregator (DTSP-версия DSUFilter)
size_t RegionAggregator::findParent(size_t x) {
    parent_[x] = (parent_[x] == x ? x : findParent(parent_[x]));
    return parent_[x];
}

void RegionAggregator::uniteSets(size_t x, size_t y) {
    x = findParent(x);
    y = findParent(y);
    if (x == y) {
        return;
    }
    for (size_t j = 0; j < values_[x].size(); ++j) {
        values_[y][j] = Interval{
            std::min(values_[y][j].lower(), values_[x][j].lower()),
            std::max(values_[y][j].upper(), values_[x][j].upper())
        };
    }
    if (componentSize_[x] > componentSize_[y]) {
        std::swap(x, y);
    }
    parent_[x] = y;
    componentSize_[y] += componentSize_[x];
}

RegionAggregator::RegionAggregator() {}

void RegionAggregator::operator()(std::vector<RadiiRegion> &elements) {
    componentSize_.resize(elements.size());
    parent_.resize(elements.size());
    values_.resize(elements.size());

    for (size_t i = 0; i < elements.size(); ++i) {
        componentSize_[i] = 1;
        parent_[i] = i;
        values_[i] = elements[i].GetIntervals();
    }

    for (size_t i = 0; i < elements.size(); ++i) {
        auto edge = values_[i];
        for (size_t j = 0; j < edge.size(); ++j) {
            auto interval = edge[j];

            // lower endpoint as key
            edge[j] = Interval{ interval.lower(), interval.lower() };
            auto it = edgesMap_.find(edge);
            if (it != edgesMap_.end()) {
                uniteSets(i, it->second);
            } else {
                edgesMap_[edge] = i;
            }

            // upper endpoint as key
            edge[j] = Interval{ interval.upper(), interval.upper() };
            it = edgesMap_.find(edge);
            if (it != edgesMap_.end()) {
                uniteSets(i, it->second);
            } else {
                edgesMap_[edge] = i;
            }

            // restore original interval
            edge[j] = interval;
        }
    }

    elements.clear();
    for (size_t i = 0; i < parent_.size(); ++i) {
        if (parent_[i] == i) {
            elements.emplace_back(values_[i]);
        }
    }
}

/// RegionExplorer (DTSP-версия Searcher)
RegionExplorer::RegionExplorer(std::vector<RadiiRegion> &results, BaseType lowerBound, BaseType upperBound)
    : resultsRef_(results), lowerBound_(lowerBound), upperBound_(upperBound) {}

void RegionExplorer::processRegion(const RadiiRegion &region,
                                   std::vector<RadiiRegion> &outRegions,
                                   std::optional<ConnectivityGraph> &graph) {
    // Если регион не слишком широкий — подготавливаем (или уточняем) connectivity graph
    if (!region.IsTooWide(upperBound_)) {
        if (!graph.has_value()) {
            OperatorLookupTable lookup_table(region.GetIntervals());
            graph.emplace(lookup_table);

            if (graph->HasOverflow()) {
                graph.reset();
            } else {
                if (!graph->IsViable()) {
                    graph.reset();
                    return;
                }
            }
        } else {
            OperatorLookupTable lookup_table(region.GetIntervals());
            graph->Refine(lookup_table);
            if (!graph->IsViable()) {
                graph->Restore();
                return;
            }
        }
    }

    // Если регион достаточно узкий — делаем "дорогую" проверку
    if (region.IsNarrowEnough(lowerBound_)) {
        if (expensiveCheck(region.GetIntervals())) {
            outRegions.emplace_back(region.GetIntervals());
        }
        if (graph.has_value()) {
            if (!graph->Restore()) {
                graph.reset();
            }
        }
        return;
    }

    // Иначе делим регион и рекурсивно обрабатываем детей
    std::vector<RadiiRegion> childrenRegions;
    region.GridSplit(childrenRegions, 2);
    for (auto &cr : childrenRegions) {
        processRegion(cr, outRegions, graph);
    }

    if (graph.has_value()) {
        if (!graph->Restore()) {
            graph.reset();
        }
    }
}

bool RegionExplorer::expensiveCheck(const RadiiRegion & /*region*/) {
    // TODO: реализация дорогой проверки (оставлено как в оригинале)
    return true;
}

void RegionExplorer::startProcessing(const std::vector<Interval> &intervals, size_t concurrency) {
    // Сортируем интервалы для стабильности
    std::vector<Interval> sortedIntervals = intervals;
    std::sort(sortedIntervals.begin(), sortedIntervals.end(),
              [](const Interval &a, const Interval &b) { return cerlt(a, b); });

    RadiiRegion rootRegion(sortedIntervals);

    // Разбиваем корневой регион на первичные задачи (1 * concurrency)
    std::vector<RadiiRegion> primaryChildren;
    rootRegion.GridSplit(primaryChildren, 1 * concurrency);

    std::vector<std::vector<RadiiRegion>> threadResults(concurrency);
    std::vector<std::thread> threads;
    std::atomic<size_t> taskIndex{0};

    for (size_t i = 0; i < concurrency; ++i) {
        threads.emplace_back(
            [&taskIndex, &threadResults, &primaryChildren, this, i] {
                while (true) {
                    size_t index = taskIndex.fetch_add(1, std::memory_order_relaxed);
                    if (index >= primaryChildren.size()) {
                        break;
                    }

                    auto progress = ((index + 1) * 100'000) / primaryChildren.size();
                    std::cerr << "\r" + std::to_string(progress / 1000) + "." +
                                     (progress % 1000 < 10 ? "0" : "") +
                                     (progress % 1000 < 100 ? "0" : "") +
                                     std::to_string(progress % 1000) + "% ";

                    std::optional<ConnectivityGraph> graph = std::nullopt;
                    processRegion(primaryChildren[index], threadResults[i], graph);
                    std::cerr << "x";
                }
            }
        );
    }

    size_t resultsSize = 0;
    for (size_t i = 0; i < concurrency; ++i) {
        threads[i].join();
        resultsSize += threadResults[i].size();
    }

    resultsRef_.clear();
    resultsRef_.reserve(resultsSize);
    for (auto &tr : threadResults) {
        resultsRef_.insert(resultsRef_.end(), tr.begin(), tr.end());
    }
    std::cerr << "\ninitial result size:\t" << resultsRef_.size() << "\n";

    // Свести соседние маленькие регионы в большие компоненты
    RegionAggregator{}(resultsRef_);

    // Отсортировать и вернуть
    std::sort(resultsRef_.begin(), resultsRef_.end(),
              [](const RadiiRegion &a, const RadiiRegion &b) {
                  return RadiiCompare{}(a.GetIntervals(), b.GetIntervals());
              });
}

} // namespace diskpack