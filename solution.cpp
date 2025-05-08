#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

const double PI = 3.141592653589793;
const double TOLERANCE = 0.15; // допустимая погрешность при проверке углов
const double STEP = 0.005; // шаг перебора значений радиусов

// Угол между двумя внешними дисками, касающимися центрального
double angle(double R, double a, double b) {
    double ra = R + a;
    double rb = R + b;
    double lab = a + b;
    double cos_theta = (ra * ra + rb * rb - lab * lab) / (2 * ra * rb);
    if (cos_theta < -1.0) cos_theta = -1.0;
    else if (cos_theta > 1.0) cos_theta = 1.0;
    return std::acos(cos_theta);
}

// Проверка замыкания кольца: сумма углов ≈ 2π
bool is_compact_ring(double R, const std::vector<double>& ring) {
    double sum = 0.0;
    int n = ring.size();
    for (int i = 0; i < n; ++i) {
        double a = ring[i];
        double b = ring[(i + 1) % n];
        sum += angle(R, a, b);
    }
    return std::abs(sum - 2 * PI) < TOLERANCE;
}

// Проверка всех допустимых конфигураций
bool is_valid_triple(double R, double r, double s, double t) {
    std::vector<std::vector<double>> patterns = {
        {r, s, t, r, s, t},
        {r, r, s, s, t, t},
        {r, s, r, s, t, t},
        {r, s, r, t, s, t},
        {r, r, r, s, t, t},
        {s, s, s, r, r, t}
    };
    for (const auto& pattern : patterns) {
        if (is_compact_ring(R, pattern)) return true;
    }
    return false;
}

int main() {
    const double R = 1.0;

    std::ofstream fout("output.txt");
    if (!fout.is_open()) {
        std::cerr << "Ошибка открытия output.txt\n";
        return 1;
    }

    fout << std::fixed << std::setprecision(5);
    fout << "Valid (r, s, t) triples for compact packing:\n";

    int count = 0;
    for (double r = 0.99; r >= 0.9; r -= STEP) {
        for (double s = r - STEP; s >= 0.9; s -= STEP) {
            for (double t = s - STEP; t >= 0.9; t -= STEP) {
                if (is_valid_triple(R, r, s, t)) {
                    fout << "(" << r << ", " << s << ", " << t << ")\n";
                    ++count;
                }
            }
        }
    }

    fout << "\nTotal valid triples: " << count << "\n";
    fout.close();

    std::cout << "Done. Found " << count << " valid triples. See output.txt\n";
    return 0;
}
