#ifndef DISKPACK_GEOMETRY_H
#define DISKPACK_GEOMETRY_H

#include <cmath>
#include <utility>
#include <limits>
#include <iostream>

namespace diskpack {

// Базовый тип для чисел с плавающей точкой
using BaseType = long double;

// Структура для хранения точки на плоскости
struct Point {
    BaseType x;
    BaseType y;

    Point() : x(0), y(0) {}
    Point(BaseType xVal, BaseType yVal) : x(xVal), y(yVal) {}
};

// Пара интервалов [low, high]
struct Interval {
    BaseType low;
    BaseType high;

    Interval() : low(0), high(0) {}
    Interval(BaseType l, BaseType h) : low(l), high(h) {}

    // Проверка корректности интервала
    bool valid() const { return low <= high; }

    // Длина интервала
    BaseType width() const { return high - low; }

    // Проверка пересечения интервалов
    bool overlaps(const Interval& other) const {
        return !(high < other.low || other.high < low);
    }

    // Сдвиг интервала
    Interval shifted(BaseType delta) const {
        return Interval(low + delta, high + delta);
    }
};

// Интервальная пара (например, для полярных координат)
struct IntervalPair {
    Interval first;
    Interval second;

    IntervalPair() = default;
    IntervalPair(const Interval& a, const Interval& b)
        : first(a), second(b) {}
};

// Представление диска на плоскости
struct Disk {
    Point center;
    BaseType radius;
    size_t typeId;

    Disk() : center(), radius(0), typeId(0) {}
    Disk(const Point& c, BaseType r, size_t t)
        : center(c), radius(r), typeId(t) {}
};

// Геометрические утилиты
namespace geometry {

inline BaseType distance(const Point& a, const Point& b) {
    BaseType dx = a.x - b.x;
    BaseType dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

inline BaseType sqDistance(const Point& a, const Point& b) {
    BaseType dx = a.x - b.x;
    BaseType dy = a.y - b.y;
    return dx * dx + dy * dy;
}

inline bool circlesTouch(const Disk& d1, const Disk& d2, BaseType eps = 1e-12) {
    BaseType dist = distance(d1.center, d2.center);
    return std::abs(dist - (d1.radius + d2.radius)) <= eps;
}

} // namespace geometry

} // namespace diskpack

#endif // DISKPACK_GEOMETRY_H