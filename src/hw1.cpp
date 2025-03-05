#include "hw1.h"
#include <cstddef>
#include <random>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

Matrix algebra::zeros(size_t n, size_t m) {
    return {n, std::vector<double>(m, 0)};
}

Matrix algebra::ones(size_t n, size_t m) {
    return {n, std::vector<double>(m, 1)};
}

Matrix algebra::random(size_t n, size_t m, double min, double max) {
    if (min > max) {
        throw std::logic_error("min cannot be greater than max");
    }

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

Matrix algebra::multiply(const Matrix &matrix, double c) {
    Matrix result(matrix);
    for (auto &row: result) {
        for (auto &elem: row) {
            elem *= c;
        }
    }
    return result;
}


Matrix algebra::multiply(const Matrix &matrix1, const Matrix &matrix2) {
    // matrix1 m x n -> 2 x 3
    // matrix2 n x p -> 3 x 4

    if (matrix1.empty() || matrix2.empty()) {
        return {};
    }

    if (matrix1[0].size() != matrix2.size()) {
        throw std::logic_error("matrices with wrong dimensions cannot be multiplied");
    }

    size_t m = matrix1.size();
    size_t n = matrix2.size();
    size_t p = matrix2[0].size();

    Matrix result(m, std::vector<double>(p));

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < p; j++) {
            for (size_t k = 0; k < n; k++) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

Matrix algebra::sum(const Matrix &matrix, double c) {
    Matrix result(matrix);

    for (auto &row: result) {
        for (auto &elem: row) {
            elem += c;
        }
    }

    return result;
}

Matrix algebra::sum(const Matrix &matrix1, const Matrix &matrix2) {
    if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
        throw std::logic_error("matrices with wrong dimensions cannot be summed");
    }

    if (matrix1.empty() || matrix2.empty()) {
        return {};
    }

    size_t m = matrix1.size();
    size_t n = matrix1[0].size();

    Matrix result(m, std::vector<double>(n));

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }

    return result;
}

Matrix algebra::transpose(const Matrix &matrix) {
    size_t m = matrix.size();
    size_t n = matrix[0].size();

    Matrix result(n, std::vector<double>(m));

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}

Matrix algebra::minor(const Matrix &matrix, size_t n, size_t m) {
    if (matrix.size() != matrix[0].size()) {
        throw std::logic_error("non-square matrices don't have minors");
    }

    size_t rows = matrix.size();
    
    if (rows == 0 || n >= rows || m >= rows) {
        throw std::logic_error("invalid input: index out of range");
    }

    // Create a matrix with dimensions (rows-1) x (rows-1)
    Matrix result(rows-1, std::vector<double>(rows-1));
    
    // Variables to track position in the result matrix
    size_t rRow = 0;
    
    // Copy elements from original matrix to result, skipping row n and column m
    for (size_t i = 0; i < rows; i++) {
        if (i == n) continue; // Skip the specified row
        
        size_t rCol = 0;
        for (size_t j = 0; j < rows; j++) {
            if (j == m) continue; // Skip the specified column
            
            // Copy the element to the result matrix
            result[rRow][rCol] = matrix[i][j];
            rCol++;
        }
        rRow++;
    }
    
    return result;
}




















