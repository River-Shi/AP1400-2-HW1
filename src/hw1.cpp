#include "hw1.h"
#include <random>
#include <iomanip>
#include <iostream>

Matrix algebra::zeros(size_t n, size_t m) {
   return {n, std::vector<double>(m, 0)};
}

Matrix algebra::ones(size_t n, size_t m) {
    return {n, std::vector<double>(m, 1)};
}

Matrix algebra::random(size_t n, size_t m, double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);

    Matrix matrix(n, std::vector<double>(m));
    for (auto &row: matrix) {
        for (auto &elem: row) {
            elem = dis(gen);
        }
    }
    return matrix;
}

void algebra::show(const Matrix &matrix) {
    for (const auto &row: matrix) {
        for (const auto &elem: row) {
            std::cout << std::fixed << std::setprecision(3) << std::setw(9) << std::right << elem << " ";
        }
        std::cout << std::endl;
    }
}
