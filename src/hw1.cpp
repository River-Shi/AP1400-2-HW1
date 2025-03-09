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
    if (matrix1.empty() && matrix2.empty()) {
        return {};
    }


    if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
        throw std::logic_error("matrices with wrong dimensions cannot be summed");
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
    if (matrix.empty()) {
        return {};
    }

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

double algebra::determinant(const Matrix &matrix) {
    if (matrix.empty()) {
        return 1;
    }

    if (matrix.size() != matrix[0].size()) {
        throw std::logic_error("non-square matrices don't have determinants");
    }

    if (matrix.size() == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }

    double det = 0;

    for (size_t i = 0; i < matrix.size(); i++) {
        double cofactor = matrix[0][i] * ((i % 2 == 0) ? 1 : -1);
        det += cofactor * determinant(minor(matrix, 0, i));
    }

    return det;
}

Matrix algebra::inverse(const Matrix &matrix) {
    if (matrix.empty()) {
        return {};
    }

    if (matrix.size() != matrix[0].size()) {
        throw std::logic_error("non-square matrices don't have inverses");
    }

    double det = determinant(matrix);

    if (std::abs(det) < 1e-10) {
        throw std::logic_error("singular matrices don't have inverses");
    }

    size_t size = matrix.size();

    // 创建一个临时矩阵存储代数余子式
    Matrix cofactors(size, std::vector<double>(size));

    // 计算所有位置的代数余子式
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            // 获取余子式矩阵 (删除第i行第j列)
            Matrix submatrix = minor(matrix, i, j);
            
            // 计算代数余子式
            double cofactorValue = determinant(submatrix);
            if ((i + j) % 2 == 1) {  // 奇数位置取负
                cofactorValue = -cofactorValue;
            }
            
            cofactors[i][j] = cofactorValue;
        }
    }

    // 计算伴随矩阵 (代数余子式矩阵的转置)
    Matrix adj = transpose(cofactors);

    // 计算逆矩阵
    Matrix result(size, std::vector<double>(size));
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            result[i][j] = adj[i][j] / det;
        }
    }

    return result;
}

Matrix algebra::concatenate(const Matrix &matrix1, const Matrix &matrix2, int axis) {
    // axis can only be 0 or 1
    if (axis == 0) {
        size_t col1 = matrix1[0].size();
        size_t col2 = matrix2[0].size();

        if (col1 != col2) {
            throw std::logic_error("matrices with wrong dimensions cannot be concatenated");
        }


        Matrix result(matrix1.size() + matrix2.size() , std::vector<double>(col1));

        for (size_t i = 0; i < matrix1.size() + matrix2.size(); i++) {
            for (size_t j = 0; j < col1; j++) {
                if (i < matrix1.size()) {
                    result[i][j] = matrix1[i][j];
                } else {
                    result[i][j] = matrix2[i-matrix1.size()][j];
                }
            }
        }

        return result;

    } else if (axis == 1) {
        size_t row1 = matrix1.size();
        size_t row2 = matrix2.size();

        if (row1 != row2) {
            throw std::logic_error("matrices with wrong dimensions cannot be concatenated");
        }

        Matrix result(row1, std::vector<double>(matrix1[0].size() + matrix2[0].size()));

        for (size_t i = 0; i < row1; i++) {
            for (size_t j = 0; j < matrix1[0].size() + matrix2[0].size(); j++) {
                if (j < matrix1[0].size()) {
                    result[i][j] = matrix1[i][j];
                } else {
                    result[i][j] = matrix2[i][j-matrix1[0].size()];
                }
            }
        }

        return result;

    } else {
        throw std::logic_error("axis can only be 0 or 1");
    }


}

Matrix algebra::ero_swap(const Matrix &matrix, size_t r1, size_t r2) {
    size_t rows = matrix.size();

    if (rows == 0 || r1 >= rows || r2 >= rows) {
        throw std::logic_error("r1 or r2 inputs are out of range");
    }

    Matrix result(matrix);
    if (r1 != r2) {
        std::swap(result[r1], result[r2]);
    }

    return result;
}

Matrix algebra::ero_multiply(const Matrix &matrix, size_t r, double c) {
    if (matrix.size() == 0 || r >= matrix.size()) {
        throw std::logic_error("r input is out of range");
    }

    Matrix result(matrix);
    for (size_t i = 0; i < matrix[r].size(); i++) {
        result[r][i] *= c;
    }

    return result;
}

Matrix algebra::ero_sum(const Matrix &matrix, size_t r1, double c, size_t r2) {
    if (matrix.size() == 0 || r1 >= matrix.size() || r2 >= matrix.size()) {
        throw std::logic_error("r1 or r2 inputs are out of range");
    }

    Matrix result(matrix);

    for (size_t i = 0; i < matrix[r1].size(); i++) {
        result[r2][i] += c * matrix[r1][i];
    }

    return result;
}

Matrix algebra::upper_triangular(const Matrix& matrix) {
    if (matrix.empty()) {
        return {};
    }

    if (matrix.size() != matrix[0].size()) {
        throw std::logic_error("non-square matrices don't have upper triangular forms");
    }

    Matrix result(matrix);
    for (size_t i = 0; i < matrix.size() - 1; i++) {
        if (std::abs(result[i][i]) < 1e-10) {
            bool find = false;
            // find the first row with non-zero value in the i-th column

            for (size_t j = i + 1; j < matrix.size(); j++) {
                if (std::abs(result[j][i]) > 1e-10) {
                    result = ero_swap(result, i, j);
                    find = true;
                    break;
                }
            }

            if (!find) {
                continue;
            }
        }

        // elimination
        for (size_t j = i + 1; j < matrix.size(); j++) {
            double c = result[j][i] / result[i][i];
            result = ero_sum(result, i, -c, j);
        }
    }

    return result;
}

















