#include <iomanip>
#include <kmdiff/linear_model.hpp>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

namespace kmdiff {

  size_t nrows(const matrix_t& m)
  {
    return m.size();
  }

  size_t ncols(const matrix_t& m)
  {
    return m[0].size();
  }

  void print_matrix(const matrix_t& m)
  {
    for (auto& i: m)
    {
      for (auto& j: i)
      {
        std::cerr << j << "\t";
      }
      std::cerr << "\n";
    }
  }

  std::string str_matrix(const matrix_t& m)
  {
    std::stringstream ss;
    ss.precision(2);
    ss << std::scientific;

    for (auto& i: m)
    {
      for (auto& j: i)
      {
        ss << j << "\t";
      }
      ss << "\n";
    }
    return ss.str();
  }

  bool is_equal_v(const vector_t& v1, const vector_t& v2)
  {
    return std::equal(v1.begin(), v1.end(), v2.begin());
  }

  bool is_equal_m(const matrix_t& m1, const matrix_t& m2)
  {
    for (size_t i=0; i<nrows(m1); i++)
      if (!is_equal_v(m1[i], m2[i]))
        return false;
    return true;
  }

  bool is_equal_d(double x, double y, double e)
  {
    return std::fabs(x-y) < e;
  }

  matrix_t transpose(const matrix_t& m)
  {
    matrix_t trp (ncols(m), vector_t(nrows(m)));
    for (size_t i=0; i<nrows(m); i++)
      for (size_t j=0; j<ncols(m); j++)
        trp[j][i] = m[i][j];
    return trp;
  }

  matrix_t multiply(const matrix_t& m1, const matrix_t& m2)
  {
    matrix_t res (nrows(m1), vector_t(ncols(m2), 0.0));
    assert(ncols(m1) == nrows(m2));
    for (size_t i=0; i<nrows(m1); i++)
      for (size_t j=0; j<ncols(m2); j++)
        for (size_t k=0; k<ncols(m1); k++)
          res[i][j] = res[i][j] + m1[i][k] * m2[k][j];
    return res;
  }

  /*
    The following code is adapted from https://github.com/atifrahman/HAWK
    https://doi.org/10.7554/eLife.32920.001
    https://doi.org/10.1371/journal.pone.0245058
  */

  std::tuple<matrix_t, matrix_t> lu_decomposition(const matrix_t& orig,
                                                  size_t n)
  {
    matrix_t lower(n, vector_t(n, 0));
    matrix_t upper(n, vector_t(n, 0));

    for (size_t i=0; i<n; i++)
    {
      // upper triangular
      for (size_t k=i; k<n; k++)
      {
        double sum = 0.0;
        for (size_t j=0; j<i; j++)
        {
          sum += (lower[i][j] * upper[j][k]);
        }
        upper[i][k] = orig[i][k] - sum;
      }
      // lower triangular
      for (size_t k=i; k<n; k++)
      {
        if (i == k)
        {
          lower[i][i] = 1;
        }
        else
        {
          double sum = 0;
          for (int j=0; j<i; j++)
          {
            sum += (lower[k][j] * upper[j][i]);
          }
          lower[k][i] = (orig[k][i] - sum) / upper[i][i];
        }
      }

    }
    return std::make_tuple(lower, upper);
  }

  std::tuple<matrix_t, bool, bool> inverse(const matrix_t& m, size_t n)
  {
    auto [lower, upper] = lu_decomposition(m, n);
    matrix_t inv(n, vector_t(n, 0));
    double det = 1;

    for (size_t invc=0; invc<n; invc++)
    {
      vector_t b(n, 0);
      for (size_t j=0; j<n; j++)
      {
        if (invc == j) b[j] = 1;
        else b[j] = 0;
      }
      vector_t y(n, 0);
      det *= lower[0][0];
      y[0] = b[0];
      for (size_t row=1; row<n; row++)
      {
        double sum = 0;
        for (size_t col=0; col < n; col++)
        {
          sum += lower[row][col] * y[col];
        }
        y[row] = b[row] - sum;
        det *= lower[row][row];
      }
      vector_t x(n, 0);
      x[n-1] = y[n-1] / upper[n-1][n-1];
      det *= upper[n-1][n-1];
      for (int row = n-2; row>-1; row--)
      {
        double sum = 0;
        for (int col = row+1; col<n; col++)
        {
          sum += upper[row][col] * x[col];
        }
        x[row] = (y[row] - sum) / upper[row][row];

        det *= upper[row][row];
      }

      for (size_t j=0; j<n; j++)
      {
        inv[j][invc] = x[j];
      }
    }

    bool singular = false, nan = false;
    if (det == 0)
      singular = true;
    else if (std::isnan(det))
      nan = true;

    return std::make_tuple(inv, singular, nan);
  }

  double sigmoid(double x)
  {
    double e = M_E;
    return 1.0 / (1.0 + std::pow(e, -x));
  }

  double linear_predictor(const vector_t& model, const vector_t& data)
  {
    double s = 0.0;
    for (size_t i=0; i<model.size(); i++)
      s += model[i]*data[i];
    return s;
  }

  double predict(const vector_t& model, const vector_t& data)
  {
    double s = 0.0;
    for (size_t i=0; i<model.size(); i++)
      s += model[i]*data[i];
    return sigmoid(s);
  }

  std::tuple<vector_t, bool, bool, double, int> glm_newton_raphson(const matrix_t& x,
                                                                   const vector_t& y,
                                                                   double gamma,
                                                                   int max_iters)
  {
    bool ise = false, ine = false;
    double re = 0.0;
    int ein = 0;
    double epsilon = 1e-6;
    double error = 0.0;
    int iter = 0;
    matrix_t A = x;
    vector_t weight_old(ncols(A), 0);
    for (size_t i=0; i<weight_old.size(); i++)
    {
      double mxx = -10000000000.0;
      for (size_t j=0; j<nrows(A); j++)
      {
        mxx = std::max(mxx, A[j][i]);
      }
      weight_old[i] = 1.0/mxx;
    }
    double prev_error = 1e18;
    matrix_t AT = transpose(x);
    matrix_t A_minus_Y(nrows(A), vector_t(1, 0));
    vector_t B_proxy(nrows(A), 0);
    matrix_t ATB(nrows(AT), vector_t(B_proxy.size(), 0));

    while (true)
    {
      double error = 0.0;
      for (size_t i=0; i<nrows(A); i++)
      {
        double z_i = 0.0;
        for (size_t j=0; j<ncols(A); j++)
        {
          z_i += A[i][j] * weight_old[j];
        }
        double alph_i = sigmoid(z_i);
        error += (y[i]-alph_i)*(y[i]-alph_i);
        B_proxy[i] = alph_i*(1.0-alph_i);
        A_minus_Y[i][0] = alph_i-y[i];
      }
      error /= nrows(A);
      re = error;
      if (std::fabs(error-prev_error) < epsilon)
        break;

      prev_error = error;
      for (size_t i=0; i<nrows(ATB); i++)
      {
        for (size_t j=0; j<ncols(ATB); j++)
        {
          ATB[i][j] = AT[i][j]*B_proxy[j];
        }
      }

      matrix_t toinv = multiply(ATB, A);

      auto [hinv, sing, nan] = inverse(toinv, nrows(toinv));
      if (sing || nan)
      {
        return make_tuple(weight_old, ise, ine, re, iter);
      }

      matrix_t gradient = multiply(AT, A_minus_Y);
      matrix_t grad_mul_hinv = multiply(hinv, gradient);

      for (size_t j=0; j<weight_old.size(); j++)
      {
        weight_old[j] -= gamma*grad_mul_hinv[j][0];
      }

      iter += 1;
      ein = iter;

      if (iter >= max_iters)
        break;

      prev_error = error;
      re = prev_error;
    }
    return make_tuple(weight_old, ise, ine, re, iter);
  }

  std::tuple<vector_t, bool, bool, double, int> glm_irls(const matrix_t& x,
                                                         const vector_t& y,
                                                         int max_iters)
  {
    double epsilon = 1e-6;
    int iter = 0;
    int ein = 0;
    bool ise = false, ine = false;
    double ret_error;
    vector_t weight(ncols(x), 1);
    matrix_t w(ncols(x), vector_t(1, 1));
    matrix_t X(nrows(x), vector_t(ncols(x), 0));
    matrix_t eta(nrows(x), vector_t(1, 0));
    matrix_t mu(nrows(x), vector_t(1, 0));
    matrix_t S(nrows(x), vector_t(1, 0));
    matrix_t z(nrows(x), vector_t(1, 0));
    matrix_t S_X(nrows(x), vector_t(ncols(x), 0));
    matrix_t S_z(nrows(x), vector_t(1, 0));

    for (size_t i=0; i<nrows(mu); i++)
    {
      mu[i][0] = (y[i] + 0.5) / 2;
      eta[i][0] = log(mu[i][0]/(1-mu[i][0]));
    }
    double prev_error = 1e18;

    while (true)
    {
      X.clear();
      S.clear();
      S_X.clear();
      z.clear();
      S_z.clear();

      double error = 0.0;
      int num_of_good = 0;
      for (size_t i=0; i<nrows(x); i++)
      {
        double g_i = mu[i][0] * (1.0 - mu[i][0]);
        if (g_i > 1e-305)
        {
          num_of_good++;
          X.push_back(x[i]);
          z.push_back(vector_t(1, eta[i][0] + (y[i] - mu[i][0]) / (g_i + 1e-305)));
          S.push_back(vector_t(1, g_i));
        }
        error += (y[i] - mu[i][0]) * (y[i] - mu[i][0]);
      }
      if (num_of_good == 0)
        break;

      error /= nrows(x);
      ret_error = error;

      if (std::fabs(error - prev_error) < epsilon)
        break;

      prev_error = error;
      matrix_t X_T = transpose(X);
      for (size_t i=0; i<nrows(S); i++)
      {
        S_X.push_back(vector_t(ncols(X), 0));
        for (size_t j=0; j<ncols(X); j++)
        {
          S_X[i][j] = S[i][0] * X[i][j];
        }
      }

      matrix_t hessian = multiply(X_T, S_X);

      auto [hessian_inv, sing, nan] = inverse(hessian, nrows(hessian));
      if (sing || nan)
      {
        ise = sing;
        ine = nan;
        ret_error = prev_error;
        break;
      }

      for (size_t i=0; i<nrows(z); i++)
      {
        S_z.push_back(vector_t(1, S[i][0] * z[i][0]));
      }
      matrix_t XTSz = multiply(X_T, S_z);
      w = multiply(hessian_inv, XTSz);

      iter += 1;
      ein = iter;

      if (iter >= max_iters)
      {
        break;
      }

      prev_error = error;
      ret_error = prev_error;

      for (size_t i=0; i<weight.size(); i++)
        weight[i] = w[i][0];

      spdlog::debug("Work eta {} {} x {} {} mu {} {} w {} {}", nrows(eta), ncols(eta), nrows(x), ncols(x), nrows(mu), ncols(mu), nrows(w), ncols(w));
      for (size_t i=0; i<nrows(x); i++)
      {
        eta[i][0] = 0;
        for (size_t j=0; j<ncols(x); j++)
        {
          eta[i][0] += x[i][j] * w[j][0];
        }
        mu[i][0] = sigmoid(eta[i][0]);
      }
      spdlog::debug("crash");
    }
    return std::make_tuple(weight, ise, ine, ret_error, ein);
  }

} // end of namespace kmdiff
