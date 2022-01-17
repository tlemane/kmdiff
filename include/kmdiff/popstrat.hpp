#pragma once

// std
#include <string>
#include <fstream>
#include <map>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <mutex>
#include <random>

// ext
#include <fmt/format.h>
#include <spdlog/spdlog.h>

// int
#include <kmdiff/utils.hpp>
#include <kmdiff/kmtricks_utils.hpp>
#include <kmdiff/kmer.hpp>
#include <kmdiff/model.hpp>
#include <kmdiff/linear_model.hpp>
#include <kmdiff/spinlock.hpp>

const std::map<std::string, std::string> parfile_map = {
  {"genotypename", "gwas_eigenstratX.geno"},
  {"snpname", "gwas_eigenstratX.snp"},
  {"indivname", "gwas_eigenstratX.ind"},
  {"evecoutname", "gwas_eigenstrat.evec"},
  {"evaloutname", "gwas_eigenstrat.eval"},
  {"usenorm", "YES"},
  {"numoutlieriter", "0"},
  {"numoutevec", "10"},
};

namespace kmdiff {

void write_parfile(const std::string& path);

void write_gwas_info(const std::string& kmfof, const std::string& path,
                     size_t nb_controls, size_t nb_cases);

void write_gwas_eigenstrat_total(const std::string& path,
                                 const std::vector<std::uint64_t>& c1,
                                 const std::vector<std::uint64_t>& c2);

void run_eigenstrat_smartpca(const std::string& popstrat_dir,
                             const std::string& parfile,
                             const std::string& log, bool is_diploid);
template<size_t MAX_C>
class EigGenoFile
{
  using ctype = typename km::selectC<MAX_C>::type;
public:
  EigGenoFile(const std::string& path)
    : m_path(path), m_out(m_path, std::ios::out)
  {
    check_fstream_good(m_path, m_out);
  }
  void push(const Range<ctype>& controls, const Range<ctype>& cases)
  {
    for (auto& c : controls)
      m_out << (c > 0 ? "1\t" : "0\t");
    for (auto& c : cases)
      m_out << (c > 0 ? "1\t" : "0\t");
    m_out << "\n";
  }

  void close() { m_out.close(); }
private:
  std::string m_path;
  std::ofstream m_out;
};

template<size_t MAX_C>
using eig_geno_t = std::shared_ptr<EigGenoFile<MAX_C>>;

class EigSnpFile
{
public:
  EigSnpFile(const std::string& path)
    : m_path(path), m_out(path, std::ios::out)
  {
    check_fstream_good(m_path, m_out);
  }

  void push()
  {
    m_out << m_i << "\t1\t0.0\t0\n";
    m_i++;
  }

  void close() { m_out.close(); }

private:
  std::string m_path;
  std::ofstream m_out;
  size_t m_i {0};
};

using eig_snp_t = std::shared_ptr<EigSnpFile>;

template<size_t MAX_C>
class Sampler
{
  using count_type = typename km::selectC<MAX_C>::type;

  public:
    Sampler(eig_geno_t<MAX_C> geno, eig_snp_t snp, double v)
      : m_geno(geno), m_snp(snp), m_v(v) {}

    bool sample()
    {
      return m_dist(m_gen) < m_v;
    }

    void sample(const Range<count_type>& r1, const Range<count_type>& r2)
    {
      std::unique_lock<spinlock> lock(m_lock);

      if (!sample())
        return;

      m_geno->push(r1, r2);
      m_snp->push();
    }

  private:
    std::default_random_engine m_gen;
    std::uniform_real_distribution<double> m_dist {0.0, 1.0};

    eig_geno_t<MAX_C> m_geno {nullptr};
    eig_snp_t m_snp {nullptr};
    double m_v {0.0};

    spinlock m_lock;
};

class PopStratCorrector
{
private:
  inline static size_t s_max_iter = 100;
  inline static size_t s_pca_count = 10;
  inline static double s_learn_rate = 0.1;
  inline static double s_epsilon = 1e-30;
  inline static bool s_stand = true;
  inline static const std::string s_m = "M";
  inline static const std::string s_f = "F";
  inline static const std::string s_u = "U";

public:
#ifndef KMDIFF_DEV_MODE
  PopStratCorrector(size_t nb_controls,
                    size_t nb_cases,
                    std::vector<size_t> control_totals,
                    std::vector<size_t> case_totals,
                    size_t npc);
#else
  PopStratCorrector(size_t nb_controls,
                    size_t nb_cases,
                    std::vector<size_t> control_totals,
                    std::vector<size_t> case_totals,
                    size_t npc,
                    bool irls);
#endif
  void load_Z(const std::string& path);
  void load_Y(const std::string& path);
  void load_C(const std::string& path);
  void load_ginfo(const std::string& path);
  void init_global_features();

  template<size_t MAX_K>
  double apply(KmerSign<MAX_K>& ks)
  {
    matrix_t local_features(m_alt_global_features);

    for (size_t i=0; i<m_size; i++)
      local_features[i][m_alt_feature_count - 1] = ks.m_counts_ratio[i]/m_totals[i];

  #ifdef KMDIFF_DEV_MODE
    vector_t model;
    bool singular, nan;
    double error
    int iter;

    if (m_use_irls)
    {
      std::tie(model, singular, nan, error, iter) = glm_irls(local_features,
                                                          m_Y,
                                                          s_learn_rate,
                                                          s_max_iter);
      if (singular || nan) crash_on_error("alt model", singular, nan, error, iter);
    }
    else
    {
      std::tie(model, singular, nan, error, iter) = glm_newton_raphson(local_features,
                                                                    m_Y,
                                                                    s_learn_rate,
                                                                    s_max_iter);
      if (singular || nan) crash_on_error("alt model", singular, nan, error, iter);
    }
  #else
    auto [model, singular, nan, error, iter] = glm_newton_raphson(local_features,
                                                                  m_Y,
                                                                  s_learn_rate,
                                                                  s_max_iter);
    if (singular || nan) crash_on_error("alt model", singular, nan, error, iter);

  #endif
    double alt_likelihood = 1.0;

    for (size_t f=0; f<nrows(local_features); f++)
    {
      vector_t data(ncols(local_features));
      for (size_t i=0; i<ncols(local_features); i++)
      {
        data[i] = local_features[f][i];
      }
      double p = predict(model, data);

      if (m_Y[f] == 1)
      {
        alt_likelihood = alt_likelihood * p;
      }
      else
      {
        alt_likelihood *= (1.0 - p);
      }
    }

    double null_likelihood = 1.0;
    for (size_t f=0; f<nrows(m_null_global_features); f++)
    {
      vector_t data(ncols(m_null_global_features));
      for (size_t i=0; i<ncols(m_null_global_features); i++)
      {
        data[i] = m_null_global_features[f][i];
      }
      double p = predict(m_null_model, data);
      if (m_Y[f] == 1)
      {
        null_likelihood *= p;
      }
      else
      {
        null_likelihood *= (1.0 - p);
      }
    }

    if (null_likelihood == 0 && alt_likelihood == 0.0)
    {
      null_likelihood = 0.001;
      alt_likelihood = 1.0;
    }

    double likelihood_ratio = null_likelihood / alt_likelihood;
    double log_likelihood_ratio = -2.0 * (log(likelihood_ratio));

    if (std::isnan(alt_likelihood))
      assert(0);

    if (std::fabs(log_likelihood_ratio) < s_epsilon ||
                  log_likelihood_ratio < 0.0 ||
                  std::isnan(alt_likelihood))
    {
      log_likelihood_ratio = 0.0;
    }
    ks.m_corrected = alglib::chisquarecdistribution(1, log_likelihood_ratio);
    return ks.m_corrected;
  }

private:
  void crash_on_error(std::string model, bool is_singular, bool is_nan, double error, int it);
  void standardize();

private:
  matrix_t m_Z;
  matrix_t m_C;
  vector_t m_Y;

  size_t m_nb_controls;
  size_t m_nb_cases;
  size_t m_size {m_nb_controls + m_nb_cases};
  size_t m_npc;

  size_t m_null_feature_count {0};
  size_t m_alt_feature_count {0};
  size_t m_cov_count {0};

  std::vector<vector_t> m_null_global_features;
  std::vector<vector_t> m_alt_global_features;

  vector_t m_null_model;

  std::vector<size_t> m_totals;
  std::vector<size_t> m_case_totals;
  std::vector<size_t> m_control_totals;

  std::vector<size_t> m_control_idx;
  std::vector<size_t> m_case_idx;

  std::vector<int> m_ginfo;

  std::mutex m_mutex;

  int m_unknown_gender {static_cast<int>(m_size)};

#ifdef DEV_MODE
  bool m_use_irls;
#endif
};

using pop_corrector_t = std::shared_ptr<PopStratCorrector>;

}; // end of namespace kmdif
