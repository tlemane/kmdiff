#include <kmdiff/popstrat.hpp>
#include <kmdiff/utils.hpp>

#define KMTRICKS_PUBLIC
#include <kmtricks/io/fof.hpp>

namespace kmdiff {

void write_parfile(const std::string& path)
{
  std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
  std::string parent = fs::path(path).parent_path().string();
  for (auto& [k, v] : parfile_map)
    out << k << ": " << fmt::format(v, parent) << "\n";
}

void write_gwas_info(const std::string& kmfof, const std::string& path,
                     size_t nb_controls, size_t nb_cases)
{
  km::Fof fof(kmfof);
  std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
  std::string parent = fs::path(path).parent_path().string();
  std::string control_ind = fmt::format("{}/control.ind", parent);
  std::ofstream coind(control_ind, std::ios::out); check_fstream_good(control_ind, coind);
  std::string case_ind = fmt::format("{}/case.ind", parent);
  std::ofstream caind(case_ind, std::ios::out); check_fstream_good(control_ind, coind);

  int i = 0;
  for (auto& [id, paths, min] : fof)
  {
    if (i < nb_controls)
    {
      out << id << "\tU\t" << "Control" << "\n";
      coind << id << "\tU\t" << "Control" << "\n";
    }
    else
    {
      out << id << "\tU\t" << "Case" << "\n";
      caind << id << "\tU\t" << "Case" << "\n";
    }
    i++;
  }
}

void write_gwas_eigenstrat_total(const std::string& path,
                                 const std::vector<std::uint64_t>& c1,
                                 const std::vector<std::uint64_t>& c2)
{
  std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
  for (auto& e : c1)
    out << e << "\n";
  for (auto& e : c2)
    out << e << "\n";
}

void run_eigenstrat_smartpca(const std::string& popstrat_dir,
                             const std::string& parfile,
                             const std::string& log, bool is_diploid)
{
  std::string cpath = fs::current_path();
  fs::current_path(popstrat_dir);
  std::string smartpca_bin = command_exists(get_binary_dir(), "smartpca");
  if (is_diploid)
    exec_external_cmd(smartpca_bin, fmt::format("-p {}", parfile), log);
  else
    exec_external_cmd(smartpca_bin, fmt::format("-V p {}", parfile), log);

  if (!fs::exists(fmt::format(parfile_map.at("evecoutname"), popstrat_dir)))
    throw EigenStratError("eigenstrat/smartpca failed.");

  std::string evec_bin = command_exists(get_binary_dir(), "evec2pca.perl");
  std::string pca_file = "gwas_eigenstrat.pca";
  std::string pcs_file = "pcs.evec";
  exec_external_cmd(evec_bin, fmt::format("{} {} {} {}",
                                          10,
                                          parfile_map.at("evecoutname"),
                                          parfile_map.at("indivname"),
                                          "gwas_eigenstrat.pca").c_str());

  std::ifstream in(pca_file, std::ios::in); check_fstream_good(pca_file, in);
  std::ofstream out(pcs_file, std::ios::out); check_fstream_good(pcs_file, in);
  std::string line; std::getline(in, line);
  int skip = std::stoi(line);
  for (int i=0; i<skip; i++)
    std::getline(in, line);

  while (std::getline(in, line))
    out << line << "\n";

  fs::current_path(cpath);
}

#ifdef KMDIFF_DEV_MODE
PopStratCorrector::PopStratCorrector(size_t nb_controls,
                                     size_t nb_cases,
                                     std::vector<size_t> control_totals,
                                     std::vector<size_t> case_totals,
                                     size_t npc,
                                     bool irls)
  : m_nb_controls(nb_controls),
    m_nb_cases(nb_cases),
    m_control_totals(control_totals),
    m_case_totals(case_totals),
    m_npc(npc),
    m_use_irls(irls)
{
  std::copy(m_control_totals.begin(), m_control_totals.end(), std::back_inserter(m_totals));
  std::copy(m_case_totals.begin(), m_case_totals.end(), std::back_inserter(m_totals));
  size_t i = 0;
  for (size_t j=0; j<m_nb_controls; j++)
    m_control_idx.push_back(i); i++;
  for (size_t j=0; j<m_nb_controls; j++)
    m_case_idx.push_back(i); i++;
}
#else
PopStratCorrector::PopStratCorrector(size_t nb_controls,
                                     size_t nb_cases,
                                     std::vector<size_t> control_totals,
                                     std::vector<size_t> case_totals,
                                     size_t npc)
  : m_nb_controls(nb_controls),
    m_nb_cases(nb_cases),
    m_control_totals(control_totals),
    m_case_totals(case_totals),
    m_npc(npc)
{
  std::copy(m_control_totals.begin(), m_control_totals.end(), std::back_inserter(m_totals));
  std::copy(m_case_totals.begin(), m_case_totals.end(), std::back_inserter(m_totals));
  size_t i = 0;
  for (size_t j=0; j<m_nb_controls; j++, i++)
    m_control_idx.push_back(i);
  for (size_t j=0; j<m_nb_cases; j++, i++)
    m_case_idx.push_back(i);
}
#endif


void PopStratCorrector::crash_on_error(std::string model, bool is_singular, bool is_nan,
                                       double e, int it)
{
  if (is_singular || is_nan)
  {
    throw SingularError(fmt::format("Matrix is non-invertible ({}), s={}, n={}, e={}, i={}",
                                    model, is_singular, is_nan, e, it));
  }
}
void PopStratCorrector::load_Z(const std::string& path)
{
  std::ifstream zin(path, std::ios::in); check_fstream_good(path, zin);
  m_Z.resize(m_size, vector_t(s_pca_count, 0));
  for (size_t i=0; i<m_Z.size(); i++)
    for(size_t j=0; j<s_pca_count; j++)
      zin >> m_Z[i][j];
  spdlog::debug("\nZ Matrix: \n{}", str_matrix(m_Z));
}

void PopStratCorrector::load_Y(const std::string& path)
{
  std::ifstream indin(path, std::ios::in); check_fstream_good(path, indin);
  for (std::string line; std::getline(indin, line);)
    m_Y.push_back(bc::utils::split(line, '\t')[2] == "Case" ? 0.0 : 1.0);
  spdlog::debug("\nY Vector: \n{}", str_vector(m_Y));
}

void PopStratCorrector::load_C(const std::string& path)
{
  if (!path.empty())
  {
    std::ifstream covin(path, std::ios::in); check_fstream_good(path, covin);
    vector_t raw;
    double val;
    while (!covin.eof())
    {
      covin >> val;
      raw.push_back(val);
    }
    matrix_t ctmp(m_size, vector_t(raw.size()/m_size, 0));
    int k = 0;
    for (size_t i=0; i<m_size; i++)
    {
      for (size_t j=0; j<(raw.size()/m_size); j++)
      {
        ctmp[i][j] = raw[k];
        k++;
      }
    }
    m_cov_count = raw.size()/m_size;

    int cidx = 0;
    m_C.resize(m_size, vector_t(raw.size()/m_size, 0));
    for (size_t i=0; i<m_nb_controls; i++)
    {
      size_t idx = 0;
      for (size_t j=0; j<ctmp[m_control_idx[i]].size(); j++)
        m_C[cidx][j] = ctmp[idx][j];
      cidx++;
    }
    for (size_t i=0; i<m_nb_cases; i++)
    {
      size_t idx = 0;
      for (size_t j=0; j<ctmp[m_case_idx[i]].size(); j++)
        m_C[cidx][j] = ctmp[idx][j];
      cidx++;
    }
    spdlog::debug("\nC Matrix: \n{}", str_matrix(m_C));
  }
}

void PopStratCorrector::load_ginfo(const std::string& path)
{
  if (!path.empty())
  {
    std::ifstream gin(path, std::ios::in); check_fstream_good(path, gin);
    size_t i=0;
    std::vector<size_t> gtmp(m_size);
    for (std::string line; std::getline(gin, line);)
    {
      if (bc::utils::split(line, '\t')[2] == s_m)
      {
        gtmp[i] = 1;
        m_unknown_gender--;
      }
      else if (bc::utils::split(line, '\t')[2] == s_f)
      {
        gtmp[i] = 0;
        m_unknown_gender--;
      }
    }

    int gidx = 0;
    m_ginfo.resize(m_size, 0);
    for (size_t i=0; i<m_nb_controls; i++)
    {
      m_ginfo[gidx] = gtmp[m_control_idx[i]];
      gidx++;
    }
    for (size_t i=0; i<m_nb_cases; i++)
    {
      m_ginfo[gidx] = gtmp[m_case_idx[i]];
      gidx++;
    }
    spdlog::debug("\nG Vector: \n{}", str_vector(m_ginfo));
  }
}

void PopStratCorrector::init_global_features()
{
  m_null_feature_count = 1 + m_npc + m_cov_count + 1;
  m_alt_feature_count = 1 + m_null_feature_count;

  m_null_global_features.resize(m_size, vector_t(m_null_feature_count, 0));
  m_alt_global_features.resize(m_size, vector_t(m_alt_feature_count, 0));

  for (size_t i=0; i < m_size; i++)
  {
    m_null_global_features[i][0] = 1;
    m_alt_global_features[i][0] = 1;

    for (size_t z=0; z<m_npc; z++)
    {
      m_null_global_features[i][z+1] = m_Z[i][z];
      m_alt_global_features[i][z+1] = m_Z[i][z];
    }

    for (size_t c=0; c<m_cov_count; c++)
    {
      m_null_global_features[i][1+m_npc+c] = m_C[i][c];
      m_alt_global_features[i][1+m_npc+c] = m_C[i][c];
    }

    if (m_unknown_gender == 0)
    {
      m_null_global_features[i][1 + m_npc + m_cov_count] = m_ginfo[i];
      m_alt_global_features[i][1 + m_npc + m_cov_count] = m_ginfo[i];
      m_null_global_features[i][1 + m_npc + m_cov_count + 1] = m_totals[i];
      m_alt_global_features[i][1 + m_npc + m_cov_count + 1] = m_totals[i];
    }
    else
    {
      m_null_global_features[i][1 + m_npc + m_cov_count] = m_totals[i];
      m_alt_global_features[i][1 + m_npc + m_cov_count] = m_totals[i];
    }
  }

  if (s_stand)
    standardize();

  spdlog::debug("\nNULL GLOBAL: \n{}", str_matrix(m_null_global_features));
  spdlog::debug("\nALT GLOBAL: \n{}", str_matrix(m_alt_global_features));

#ifdef KMDIFF_DEV_MODE
  if (m_use_irls)
  {
    auto [model, singular, nan, error, iter] = glm_irls(m_null_global_features,
                                                        m_Y,
                                                        s_learn_rate,
                                                        s_max_iter);
    if (singular || nan) crash_on_error("null model", model, singular, nan, error, iter);
    m_null_model = std::move(model);
  }
  else
  {
    auto [model, singular, nan, error, iter] = glm_newton_raphson(m_null_global_features,
                                                                  m_Y,
                                                                  s_learn_rate,
                                                                  s_max_iter);
    if (singular || nan) crash_on_error("null model", singular, nan, error, iter);
    m_null_model = std::move(model);
  }
#else
  auto [model, singular, nan, error, iter] = glm_newton_raphson(m_null_global_features,
                                                                m_Y,
                                                                s_learn_rate,
                                                                s_max_iter);
  if (singular || nan) crash_on_error("null model", singular, nan, error, iter);
  m_null_model = model;
  spdlog::debug("null model {} -> {} ", m_null_model.size(), str_vector(m_null_model));
#endif
}

void PopStratCorrector::standardize()
{
  vector_t means(ncols(m_null_global_features), 0);
  vector_t stddev(nrows(m_null_global_features), 0);

  for (size_t i=0; i<nrows(m_null_global_features); i++)
    for (size_t j=0; j<ncols(m_null_global_features); j++)
      means[j] += m_null_global_features[i][j];

  for (size_t i=1; i<ncols(m_null_global_features); i++)
    means[i] /= ncols(m_null_global_features);

  for (size_t i=0; i<nrows(m_null_global_features); i++)
    for (size_t j=1; j<ncols(m_null_global_features); j++)
      stddev[j] += std::pow(m_null_global_features[i][j] - means[j], 2);

  for (size_t i=1; i<ncols(m_null_global_features); i++)
  {
    stddev[i] /= nrows(m_null_global_features);
    stddev[i] = std::sqrt(stddev[i]);
  }

  for (size_t i=0; i<nrows(m_null_global_features); i++)
  {
    for (size_t j=1; j<ncols(m_null_global_features); j++)
    {
      if (std::fabs(stddev[i]) > 1e-305)
      {
        m_null_global_features[i][j] = (m_null_global_features[i][j] - means[j])/stddev[i];
        m_alt_global_features[i][j] = (m_alt_global_features[i][j] - means[j])/stddev[i];
      }
    }
  }
}
}; // end of namespace kmdiff
