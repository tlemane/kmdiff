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
#include <kmdiff/progress.hpp>
#include <kmdiff/threadpool.hpp>
#include <kmdiff/accumulator.hpp>

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
      Sampler(eig_geno_t<MAX_C> geno, eig_snp_t snp, double v, std::size_t seed = 0)
        : m_geno(geno), m_snp(snp), m_v(v)
      {
        if (seed)
          m_gen.seed(seed);
      }

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

  class pop_strat_corrector
  {
    public:
      inline static std::size_t s_max_iter = 100;
      inline static std::size_t s_pca_count = 10;
      inline static double s_learn_rate = 0.1;
      inline static double s_epsilon = 1e-30;
      inline static bool s_stand = true;
      inline static bool s_irls = false;

      static constexpr char s_m = 'M';
      static constexpr char s_f = 'F';
      static constexpr char s_u = 'U';

    public:
      pop_strat_corrector(std::size_t nb_controls, std::size_t nb_cases,
                          const vector_ull_t& control_totals,
                          const vector_ull_t& case_totals,
                          std::size_t npc);

      void load_Z(const std::string& path);
      void load_Y(const std::string& path);
      void load_C(const std::string& path);
      void load_ginfo(const std::string& path);
      void init_global_features();

      template<size_t KSIZE>
      void apply(std::vector<acc_t<KmerSign<KSIZE>>>& accumulators,
                 std::vector<acc_t<KmerSign<KSIZE>>>& pop_accumulators,
                 std::size_t nb_threads)
      {
        ThreadPool pool(nb_threads);

        const auto size = accumulators.size();

        std::exception_ptr ep = nullptr;

        indicators::ProgressBar* pb = nullptr;

        if ((spdlog::get_level() != spdlog::level::debug) && isatty_stderr())
        {
          pb = get_progress_bar("progress", size, 50, indicators::Color::white, false);
          pb->print_progress();
        }

        for (std::size_t p = 0; p < size; p++)
        {
          auto worker = [&acc = accumulators[p], &pacc = pop_accumulators[p], this, pb, &ep](int id) {

            try
            {
              while (auto& oks = acc->get())
              {
                spdlog::debug(oks.value().to_string());

                this->apply(oks.value());
                auto coks = oks.value();
                pacc->push(std::move(coks));
              }
            } catch (...) { ep = std::current_exception(); }

            acc->destroy();
            pacc->finish();

            if (pb)
              pb->tick();
          };

          pool.add_task(worker);
        }

        pool.join_all();
        delete pb;

        if (ep != nullptr)
          rethrow_exception(ep);
      }

    private:

      void standardize();

      template<std::size_t KSIZE>
      void apply(KmerSign<KSIZE>& ks)
      {
        matrix_t local_features(m_alt_global_features);

        for (std::size_t i = 0; i < m_size; i++)
        {
          local_features[i][m_alt_feature_count - 1] = ks.m_counts_ratio[i] / m_totals[i];
        }

        #ifdef KMD_USE_IRLS
          auto [model, singular, nan, error, iter] = glm_newton_irls(
            local_features, m_Y, s_max_iter);
        #else
          auto [model, singular, nan, error, iter] = glm_newton_raphson(
            local_features, m_Y, s_learn_rate, s_max_iter);
        #endif

        double alt_likelihood = 1.0;

        for (std::size_t f = 0; f < nrows(local_features); f++)
        {
          vector_t data(ncols(local_features));

          for (std::size_t i = 0; i < ncols(local_features); i++)
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
            alt_likelihood *= 1.0 - p;
          }
        }

        double null_likelihood = 1.0;

        for (std::size_t f = 0; f < nrows(m_null_global_features); f++)
        {
          vector_t data(ncols(m_null_global_features));

          for (std::size_t i = 0; i < ncols(m_null_global_features); i++)
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
            null_likelihood *= 1.0 - p;
          }
        }

        if (null_likelihood == 0.0 && alt_likelihood == 0.0)
        {
          null_likelihood = 0.001;
          alt_likelihood = 1.0;
        }

        double likelihood_ratio = null_likelihood / alt_likelihood;
        double log_likelihood_ratio = -2.0 * (log(likelihood_ratio));

        if (std::fabs(log_likelihood_ratio) < s_epsilon ||
            log_likelihood_ratio < 0.0 ||
            std::isnan(alt_likelihood))
        {
          log_likelihood_ratio = 0.0;
        }

        auto corrected = alglib::chisquarecdistribution(1, log_likelihood_ratio);

        spdlog::debug("{} {} {} {}", ks.to_string(), ks.m_pvalue, corrected, significance_to_char(ks.m_sign));

        ks.set_pval(corrected);
      }

    private:
      matrix_t m_Z;
      matrix_t m_C;
      vector_t m_Y;

      std::size_t m_nb_controls {0};
      std::size_t m_nb_cases {0};

      std::size_t m_size {m_nb_controls + m_nb_cases};

      std::size_t m_null_feature_count {0};
      std::size_t m_alt_feature_count {0};

      std::size_t m_cov_count {0};

      vector_t m_null_model;

      vector_ull_t m_totals;
      vector_ull_t m_control_totals;
      vector_ull_t m_case_totals;

      vector_ull_t m_control_idx;
      vector_ull_t m_case_idx;

      matrix_t m_null_global_features;
      matrix_t m_alt_global_features;

      std::vector<int> m_ginfo;

      std::size_t m_unkg {m_size};

      std::size_t m_npc {0};
  };

  using pop_strat_corrector_t = std::shared_ptr<pop_strat_corrector>;

} // end of namespace kmdif

