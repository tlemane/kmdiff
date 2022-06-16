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
                       size_t nb_controls, size_t nb_cases, const std::string& gender_file)
  {
    km::Fof fof(kmfof);

    std::map<std::string, char> gm;

    if (gender_file.size() > 0)
    {
      std::ifstream in_g(path, std::ios::in);

      check_fstream_good(path, in_g);

      for (std::string line; std::getline(in_g, line);)
      {
        auto vl = bc::utils::split(line, ' ');
        fof.get_i(vl[0]);

        char g = vl[1][0];

        if (g != 'U' || g != 'M' || g != 'F')
        {
          throw IOError(fmt::format("Unknown gender: {}", g));
        }
      }
    }

    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    std::string parent = fs::path(path).parent_path().string();
    std::string control_ind = fmt::format("{}/control.ind", parent);
    std::ofstream coind(control_ind, std::ios::out); check_fstream_good(control_ind, coind);
    std::string case_ind = fmt::format("{}/case.ind", parent);
    std::ofstream caind(case_ind, std::ios::out); check_fstream_good(control_ind, coind);

    std::size_t i = 0;
    for (auto& [id, paths, min] : fof)
    {
      if (i < nb_controls)
      {
        auto it = gm.find(id);
        if (it == gm.end())
        {
          out << id << "\tU\t" << "Control" << "\n";
          coind << id << "\tU\t" << "Control" << "\n";
        }
        else
        {
          out << id << '\t' << gm.at(id) << '\t' << "Control" << "\n";
          coind << id << '\t' << gm.at(id) <<'\t' << "Control" << "\n";
        }
      }
      else
      {
        auto it = gm.find(id);
        if (it == gm.end())
        {
          out << id << "\tU\t" << "Case" << "\n";
          caind << id << "\tU\t" << "Case" << "\n";
        }
        else
        {
          out << id << '\t' << gm.at(id) << '\t' << "Case" << "\n";
          caind << id << '\t' << gm.at(id) << '\t' << "Case" << "\n";
        }
      }
      i++;
    }
  }

  void write_gwas_eigenstrat_total(const std::string& path,
                                   const std::vector<std::size_t>& c1,
                                   const std::vector<std::size_t>& c2)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    for (auto& e : c1)
      out << e << "\n";
    for (auto& e : c2)
      out << e << "\n";
  }

  void run_eigenstrat_smartpca(const std::string& popstrat_dir,
                               const std::string& parfile,
                               const std::string& log,
                               bool is_diploid,
                               std::size_t n)
  {
    std::string cpath = fs::current_path();
    fs::current_path(popstrat_dir);
    std::string smartpca_bin = command_exists(get_binary_dir(), "smartpca");
    if (is_diploid)
      exec_external_cmd(smartpca_bin, fmt::format("-p {}", parfile), log);
    else
      exec_external_cmd(smartpca_bin, fmt::format("-V -p {}", parfile), log);

    if (!fs::exists(fmt::format(parfile_map.at("evecoutname"), popstrat_dir)))
      throw EigenStratError("eigenstrat/smartpca failed.");

    std::string evec_bin = command_exists(get_binary_dir(), "evec2pca.perl");
    std::string pca_file = "gwas_eigenstrat.pca";
    std::string pcs_file = "pcs.evec";
    exec_external_cmd(evec_bin, fmt::format("{} {} {} {}",
                                            n < 10 ? n : 10,
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

  pop_strat_corrector::pop_strat_corrector(std::size_t nb_controls, std::size_t nb_cases,
                                           const vector_ull_t& control_totals,
                                           const vector_ull_t& case_totals,
                                           std::size_t npc)
    : m_nb_controls(nb_controls), m_nb_cases(nb_cases),
      m_control_totals(control_totals), m_case_totals(case_totals),
      m_npc(npc)
  {
    std::copy(m_control_totals.begin(), m_control_totals.end(), std::back_inserter(m_totals));
    std::copy(m_case_totals.begin(), m_case_totals.end(), std::back_inserter(m_totals));

    m_control_idx.reserve(m_nb_controls); m_case_idx.reserve(m_nb_cases);

    for (std::size_t i = 0; i < m_nb_controls; i++) { m_control_idx.push_back(i); }
    for (std::size_t i = 0; i < m_nb_cases; i++) { m_case_idx.push_back(i + m_nb_controls); }
  }

  void pop_strat_corrector::load_Z(const std::string& path)
  {
    std::ifstream zin(path, std::ios::in); check_fstream_good(path, zin);
    m_Z.resize(m_size, vector_t(s_pca_count, 0));
    for (std::size_t i = 0; i < m_Z.size(); i++)
      for (std::size_t j = 0; j < s_pca_count; j++)
        zin >> m_Z[i][j];

    spdlog::debug("\nZ Matrix:\n{}", str_matrix(m_Z));
  }

  void pop_strat_corrector::load_Y(const std::string& path)
  {
    std::ifstream yin(path, std::ios::in); check_fstream_good(path, yin);
    for (std::string line; std::getline(yin, line);)
      m_Y.push_back(bc::utils::split(line, '\t')[2] == "Case" ? 0.0 : 1.0);

    spdlog::debug("\nY Vector:\n{}", str_vector(m_Y));
  }

  void pop_strat_corrector::load_C(const std::string& path)
  {
    if (!path.empty())
    {
      std::ifstream coin(path, std::ios::in); check_fstream_good(path, coin);

      vector_t raw;
      double val;

      while (!coin.eof())
      {
        coin >> val;
        raw.push_back(val);
      }

      matrix_t ctmp(m_size, vector_t(raw.size()/m_size, 0));

      int k = 0;
      for (std::size_t i = 0; i < m_size; i++)
      {
        for (std::size_t j = 0; j < (raw.size() / m_size); j++)
        {
          ctmp[i][j] = raw[k];
        }
      }

      m_cov_count = raw.size() / m_size;

      int cidx = 0;
      m_C.resize(m_size, vector_t(raw.size() / m_size, 0));

      for (std::size_t i = 0; i < m_nb_controls; i++)
      {
        std::size_t idx = 0;
        for (std::size_t j = 0; ctmp[m_control_idx[i]].size(); j++)
        {
          m_C[cidx][j] = ctmp[idx][j];
        }
        cidx++;
      }

      for (std::size_t i = 0; i < m_nb_cases; i++)
      {
        std::size_t idx = 0;
        for (std::size_t j = 0; j < ctmp[m_case_idx[i]].size(); j++)
        {
          m_C[cidx][j] = ctmp[idx][j];
        }
        cidx++;
      }
    }
    spdlog::debug("\nC Matrix:\n{}", str_matrix(m_C));
  }

  void pop_strat_corrector::load_ginfo(const std::string& path)
  {
    if (!path.empty())
    {
      std::ifstream gin(path, std::ios::in); check_fstream_good(path, gin);

      std::size_t i = 0;

      vector_ull_t gtmp(m_size);

      for (std::string line; std::getline(gin, line);)
      {
        if (bc::utils::split(line, '\t')[2].at(0) == s_m)
        {
          gtmp[i] = 1;
          m_unkg--;
        }
        else if (bc::utils::split(line, '\t')[2].at(0) == s_f)
        {
          gtmp[i] = 0;
          m_unkg--;
        }
      }

      int gidx = 0;
      m_ginfo.resize(m_size, 0);

      for (std::size_t i = 0; i < m_nb_controls; i++)
      {
        m_ginfo[gidx] = gtmp[m_control_idx[i]];
        gidx++;
      }

      for (std::size_t i = 0; i < m_nb_cases; i++)
      {
        m_ginfo[gidx] = gtmp[m_case_idx[i]];
        gidx++;
      }
    }

    spdlog::debug("\nG Vector:\n{}", str_vector(m_ginfo));
  }

  void pop_strat_corrector::init_global_features()
  {
    m_null_feature_count = 1 + m_npc + m_cov_count + 1;
    m_alt_feature_count = 1 + m_null_feature_count;

    m_null_global_features.resize(m_size, vector_t(m_null_feature_count, 0));
    m_alt_global_features.resize(m_size, vector_t(m_alt_feature_count, 0));

    for (std::size_t i = 0; i < m_size; i++)
    {
      m_null_global_features[i][0] = 1;
      m_alt_global_features[i][0] = 1;

      for (std::size_t z = 0; z < m_npc; z++)
      {
        m_null_global_features[i][z + 1] = m_Z[i][z];
        m_alt_global_features[i][z + 1] = m_Z[i][z];
      }

      for (std::size_t c = 0; c < m_cov_count; c++)
      {
        m_null_global_features[i][1 + m_npc + c] = m_C[i][c];
        m_alt_global_features[i][1 + m_npc + c] = m_C[i][c];
      }

      if (m_unkg == 0)
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

    spdlog::debug("\nNULL:\n{}", str_matrix(m_null_global_features));
    spdlog::debug("\nALT:\n{}", str_matrix(m_alt_global_features));

    #ifdef KMD_USE_IRLS
      auto [model, singular, nan, error, iter] = glm_irls(
        m_null_global_features, m_Y, s_max_iter);
    #else
      auto [model, singular, nan, error, iter] = glm_newton_raphson(
        m_null_global_features, m_Y, s_learn_rate, s_max_iter);
    #endif

    m_null_model = std::move(model);
  }

  void pop_strat_corrector::standardize()
  {
    vector_t means(ncols(m_null_global_features), 0);
    vector_t stddev(nrows(m_null_global_features), 0);

    for (std::size_t i = 0; i < nrows(m_null_global_features); i++)
    {
      for (std::size_t j = 0; j < ncols(m_null_global_features); j++)
      {
        means[j] += m_null_global_features[i][j];
      }
    }

    for (std::size_t i = 1; i < ncols(m_null_global_features); i++)
    {
      means[i] /= ncols(m_null_global_features);
    }

    for (std::size_t i = 0; i < nrows(m_null_global_features); i++)
    {
      for (std::size_t j = 1; j < ncols(m_null_global_features); j++)
      {
        stddev[j] += std::pow(m_null_global_features[i][j] - means[j], 2);
      }
    }

    for (std::size_t i = 1; i < ncols(m_null_global_features); i++)
    {
      stddev[i] /= nrows(m_null_global_features);
      stddev[i] = std::sqrt(stddev[i]);
    }

    for (std::size_t i = 0; i < nrows(m_null_global_features); i++)
    {
      for (std::size_t j = 1; j < ncols(m_null_global_features); j++)
      {
        if (std::fabs(stddev[i]) > 1e-305)
        {
          m_null_global_features[i][j] = (m_null_global_features[i][j] - means[j]) / stddev[i];
          m_alt_global_features[i][j] = (m_alt_global_features[i][j] - means[j]) / stddev[i];
        }
      }
    }
  }

} // end of namespace kmdiff

