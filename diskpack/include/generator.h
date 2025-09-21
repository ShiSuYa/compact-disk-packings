#pragma once

#include <diskpack/corona.h>
#include <set>
#include <random>

namespace diskpack {

/// Результат процесса построения упаковки
enum class GenerationResult {
    Success,          ///< Упаковка построена полностью
    Impossible,       ///< С такими радиусами упаковка невозможна
    PrecisionLimit,   ///< Нарушена верхняя граница точности
    CoronaViolation   ///< Нарушено условие короны (обычно из-за слишком высокой точности)
};

/// Очередь дисков с приоритетом — ближайшие к (0,0) идут первыми
using DiskPriorityQueue = std::multiset<DiskPointer, decltype(&LessNormCompare)>;

std::ostream& operator<<(std::ostream& out, GenerationResult result);

/// Базовый генератор упаковок
class PackingBuilderBase {
protected:
    std::mt19937 rngEngine;

    std::vector<Interval> diskRadii;          ///< Радиусы дисков
    const BaseType maxPrecisionWidth;         ///< Верхняя граница ширины интервала радиусов
    BaseType targetPackingRadius;             ///< Желаемый радиус упаковки
    size_t maxDiskCount;                      ///< Максимальное число дисков в упаковке
    std::list<DiskPointer> currentPacking;    ///< Текущее состояние упаковки
    DiskPriorityQueue pendingDisks;           ///< Очередь обрабатываемых дисков
    OperatorLookupTable operatorTable;        ///< Таблица операторов для вычислений
    std::vector<size_t> diskUsageCounter;     ///< Учёт использования дисков
    size_t ignoreRadiusLimit;                 ///< Лимит на количество игнорируемых радиусов
    BaseType achievedPackingRadius;           ///< Фактический радиус достигнутой упаковки

    /// Заполнение короны вокруг диска
    GenerationResult fillCorona(Corona& corona);

    /// Выбор диска из очереди и заполнение его короны
    GenerationResult stepForward();

    /// Перемешивание индексов для рандомизации порядка перебора
    void shuffleIndices(std::vector<size_t>& order);

    /// Добавление нового диска в упаковку
    void pushDisk(Disk&& disk, size_t typeIndex);

    /// Удаление диска (LIFO)
    void popDisk(size_t typeIndex);

    /// Запомнить радиус упаковки по самому дальнему диску
    void updateAchievedRadius(const Disk& furthest);

    /// Проверка на пересечение с уже размещёнными дисками
    bool intersectsExisting(const Disk& candidate) const;

    /// Проверка, что диск полностью помещается в область
    bool fitsInsideBounds(const Disk& disk) const;

    /// Проверка соблюдения ограничений упаковки
    bool checkConstraints() const;

    /// Проверка, что упаковка уже достаточно большая
    bool largeEnough() const;

public:
    std::optional<ConnectivityGraph> connectivityMap;

    PackingBuilderBase(const std::vector<Interval>& radii,
                       const BaseType& targetRadius,
                       const BaseType& precisionLimit,
                       const size_t& diskCountLimit,
                       const size_t& maxIgnored = 0);

    /// Полная генерация упаковки с нуля
    GenerationResult generate(const size_t& centerDiskType);

    /// Сброс состояния
    void reset();

    /// Продолжение генерации после изменения параметров
    GenerationResult resume();

    /// Задать новый желаемый радиус упаковки
    void setTargetRadius(const BaseType& newRadius);

    /// Задать новый предел на число дисков
    void setDiskCountLimit(const size_t& newLimit);

    /// Обновить набор радиусов
    void setRadii(const std::vector<Interval>& newRadii);

    /// Получить фактически достигнутый радиус упаковки
    const BaseType& getAchievedRadius();

    /// Получить текущий радиус упаковки
    const BaseType& getTargetRadius();

    /// Доступ к текущему набору дисков
    const std::list<DiskPointer>& getPacking() const;
};

} // namespace diskpack