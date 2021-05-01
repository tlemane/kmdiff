#include <kmdiff/popstrat.hpp>

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
  fof_t fof = parse_km_fof(kmfof);
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

void write_gwas_eigenstrat_total(const std::string& kmdir, const std::string& path)
{
  std::string kmtotal = fmt::format("{}/storage/kmers.total", kmdir);
  std::ifstream in(kmtotal, std::ios::in); check_fstream_good(kmtotal, in);
  std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
  for (std::string line; std::getline(in, line);)
  {
    out << bc::utils::split(line, ' ')[1] << "\n";
  }
}

void run_eigenstrat_smartpca(const std::string& popstrat_dir,
                             const std::string& parfile,
                             const std::string& log, bool is_diploid)
{
  std::string cpath = fs::current_path();
  fs::current_path(popstrat_dir);
  std::string smartpca_bin = command_exists(get_binary_dir(), "smartpca");
  if (is_diploid)
    std::system(fmt::format("{} -p {} > {}", smartpca_bin, parfile, log).c_str());
  else
    std::system(fmt::format("{} -V -p {} > {}", smartpca_bin, parfile, log).c_str());
  if (!fs::exists(fmt::format(parfile_map.at("evecoutname"), popstrat_dir)))
    throw EigenStratError("eigenstrat/smartpca failed.");

  std::string evec_bin = command_exists(get_binary_dir(), "evec2pca.perl");
  std::string pca_file = "gwas_eigenstrat.pca";
  std::string pcs_file = "pcs.evec";
  std::system(fmt::format("{} {} {} {} {}",
                          evec_bin,
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

}
#endif

void PopStratCorrector::load_Z(const std::string& path)
{
  std::ifstream zin(path, std::ios::in); check_fstream_good(path, zin);
  m_Z.resize(m_size, vector_t(s_learn_rate, 0));
  for (size_t i=0; i<m_Z.size(); i++)
    for(size_t j=0; j<s_pca_count; j++)
      zin >> m_Z[i][j];
}

void PopStratCorrector::load_Y(const std::string& path)
{
  std::ifstream indin(path, std::ios::in); check_fstream_good(path, indin);
  for (std::string line; std::getline(indin, line);)
    m_Y.push_back(bc::utils::split(line, '\t')[2] == "Case" ? 0.0 : 1.0);
}

void PopStratCorrector::load_ginfo(const std::string& path)
{
  std::ifstream gin(path, std::ios::in); check_fstream_good(path, gin);
  for (std::string line; std::getline(gin, line);)
  {

  }
}

void PopStratCorrector::init_global_features()
{
  size_t cov_count;
  m_null_feature_count = 1 + m_npc + cov_count + 1;
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

    for (size_t c=0; c<cov_count; c++)
    {
      m_null_global_features[i][1+m_npc+c] = m_C[i][c];
      m_alt_global_features[i][1+m_npc+c] = m_C[i][c];
    }

    if (m_unknown_gender == 0)
    {
      m_null_global_features[i][1 + m_npc + cov_count] = m_ginfo[i];
      m_alt_global_features[i][1 + m_npc + cov_count] = m_ginfo[i];
      m_null_global_features[i][1 + m_npc + cov_count + 1] = m_totals[i];
      m_alt_global_features[i][1 + m_npc + cov_count + 1] = m_totals[i];
    }
    else
    {
      m_null_global_features[i][1 + m_npc + cov_count] = m_totals[i];
      m_alt_global_features[i][1 + m_npc + cov_count] = m_totals[i];
    }
  }

#ifdef KMDIFF_DEV_MODE
  if (m_use_irls)
  {
    auto [model, singular, nan, error, iter] = glm_irls(m_alt_global_features,
                                                        m_Y,
                                                        s_learn_rate,
                                                        s_max_iter);
  }
  else
  {
    auto [model, singular, nan, error, iter] = glm_newton_raphson(m_alt_global_features,
                                                                  m_Y,
                                                                  s_learn_rate,
                                                                  s_max_iter);
  }
#else
  auto [model, singular, nan, error, iter] = glm_newton_raphson(m_alt_global_features,
                                                                m_Y,
                                                                s_learn_rate,
                                                                s_max_iter);

#endif

  m_null_model = std::move(model);
}

double PopStratCorrector::apply(vector_t& counts_ratio)
{
  matrix_t local_features(m_alt_global_features);

  for (size_t i=0; i<m_size; i++)
    local_features[i][m_alt_feature_count - 1] = counts_ratio[i];

#ifdef KMDIFF_DEV_MODE
  if (m_use_irls)
  {
    auto [model, singular, nan, error, iter] = glm_irls(local_features,
                                                        m_Y,
                                                        s_learn_rate,
                                                        s_max_iter);
  }
  else
  {
    auto [model, singular, nan, error, iter] = glm_newton_raphson(local_features,
                                                                  m_Y,
                                                                  s_learn_rate,
                                                                  s_max_iter);
  }
#else
  auto [model, singular, nan, error, iter] = glm_newton_raphson(local_features,
                                                                m_Y,
                                                                s_learn_rate,
                                                                s_max_iter);

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
      alt_likelihood *= p;
    }
    else
    {
      alt_likelihood *= (1.0 - p);
    }
  }

  // why here ???? TODO
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
      alt_likelihood *= p;
    }
    else
    {
      alt_likelihood *= (1.0 - p);
    }
  }
  // I think the code above is redundant, requires careful checking

  if (null_likelihood == 0 && alt_likelihood == 0.0)
  {
    null_likelihood = 0.001;
    alt_likelihood = 1.0;
  }

  double likelihood_ratio = null_likelihood / alt_likelihood;
  double log_likelihood_ratio = -2.0 * (log(likelihood_ratio));

  if (std::fabs(log_likelihood_ratio) < s_epsilon || log_likelihood_ratio < 0.0)
    log_likelihood_ratio = 0.0;

  return alglib::chisquarecdistribution(1, log_likelihood_ratio);
}

}; // end of namespace kmdiff