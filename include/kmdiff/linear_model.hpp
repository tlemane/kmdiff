#pragma once
#include <vector>
#include <tuple>

#include <xmmintrin.h>

namespace kmdiff {

using matrix_t = std::vector<std::vector<double>>;
using vector_t = std::vector<double>;

size_t nrows(const matrix_t& m);
size_t ncols(const matrix_t& m);
void print_matrix(const matrix_t& m);
bool is_equal_v(const vector_t& v1, const vector_t& v2);
bool is_equal_m(const matrix_t& m1, const matrix_t& m2);
bool is_equal_d(double x, double y, double e = 1e-15);

matrix_t transpose(const matrix_t& m);
matrix_t multiply(const matrix_t& m1, const matrix_t& m2);
std::tuple<matrix_t, matrix_t> lu_decomposition(const matrix_t& orig, size_t n);
std::tuple<matrix_t, bool, bool> inverse(const matrix_t& m, size_t n);

double sigmoid(double x);
double linear_predictor(const vector_t& model, const vector_t& data);
double predict(const vector_t& model, const vector_t& data);
std::tuple<vector_t, bool, bool, double, int> glm_newton_raphson(const matrix_t& x,
                                                                 const vector_t& y,
                                                                 double gamma,
                                                                 int max_iters);
std::tuple<vector_t, bool, bool, double, int> glm_irls(const matrix_t& x,
                                                       const vector_t& y,
                                                       double gamma,
                                                       int max_iters);
}; // end of namespace kmdiff


