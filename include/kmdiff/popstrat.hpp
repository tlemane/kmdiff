#pragma once

// std
#include <string>
#include <fstream>
#include <map>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <mutex>

// ext
#include <fmt/format.h>
#include <spdlog/spdlog.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/utilities.hpp>

// int
#include <kmdiff/utils.hpp>
#include <kmdiff/kmtricks_utils.hpp>
#include <kmdiff/kmer.hpp>
#include <kmdiff/model.hpp>
#include <kmdiff/linear_model.hpp>

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
void write_gwas_eigenstrat_total(const std::string& kmdir, const std::string& path);
void run_eigenstrat_smartpca(const std::string& popstrat_dir,
                             const std::string& parfile,
                             const std::string& log, bool is_diploid);

template<size_t MAX_C>
class EigGenoFile
{
  using ctype = typename selectC<MAX_C>::type;
public:
  EigGenoFile(const std::string& path)
    : m_path(path), m_out(m_path, std::ios::out)
  {
    check_fstream_good(m_path, m_out);
  }
  void push(const Range<ctype>& controls, const Range<ctype>& cases)
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    for (auto& c : controls)
      m_out << (c > 0 ? "1\t" : "0\t");
    for (auto& c : cases)
      m_out << (c > 0 ? "1\t" : "0\t");
    m_out << "\n";
  }

private:
  std::string m_path;
  std::ofstream m_out;
  std::mutex m_mutex;
};

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
    std::unique_lock<std::mutex> lock(m_mutex);
    m_out << std::to_string(m_i) << "\t1\t0.0\t0\n";
    m_i++;
  }

private:
  std::string m_path;
  std::ofstream m_out;
  std::mutex m_mutex;
  size_t m_i;
};

class PopStratCorrector : public ICorrector
{
private:
  inline static size_t s_max_iter = 25;
  inline static size_t s_pca_count = 10;
  inline static double s_learn_rate = 0.1;
  inline static double s_epsilon = 1e-30;

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
  void load_ginfo(const std::string& path);
  void init_global_features();

  double apply(vector_t& counts_ratio);

private:
  matrix_t m_Z;
  matrix_t m_C;
  vector_t m_Y;

  size_t m_nb_controls;
  size_t m_nb_cases;
  size_t m_size;
  size_t m_npc;

  size_t m_null_feature_count {0};
  size_t m_alt_feature_count {0};

  std::vector<vector_t> m_null_global_features;
  std::vector<vector_t> m_alt_global_features;

  vector_t m_null_model;

  std::vector<size_t> m_totals;
  std::vector<size_t> m_case_totals;
  std::vector<size_t> m_control_totals;

  std::vector<int> m_ginfo;

  int m_unknown_gender;

#ifdef DEV_MODE
  bool m_use_irls;
#endif
};

using pop_corrector_t = std::shared_ptr<PopStratCorrector>;

}; // end of namespace kmdif