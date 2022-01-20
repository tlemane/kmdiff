/*****************************************************************************
 *   kmdiff
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

// std
#include <memory>
#include <string>
#include <vector>
#include <random>

// int
#include <kmdiff/accumulator.hpp>
#include <kmdiff/model.hpp>
#include <kmdiff/threadpool.hpp>
#include <kmdiff/utils.hpp>
#include <kmdiff/blocking_queue.hpp>
#include <kmdiff/popstrat.hpp>
#include <kmdiff/range.hpp>
#include <kmdiff/time.hpp>
#include <kmdiff/progress.hpp>

#define KMTRICKS_PUBLIC
#include <kmtricks/merge.hpp>
#include <kmtricks/utils.hpp>

namespace kmdiff
{

template<std::size_t KSIZE, std::size_t CMAX>
class diff_observer : public km::IMergeObserver<KSIZE, CMAX>
{
  using count_type = typename km::selectC<CMAX>::type;
  public:
    diff_observer(const std::shared_ptr<IModel<CMAX>>& model,
                  acc_t<KmerSign<KSIZE>> acc,
                  double threshold,
                  std::size_t controls,
                  std::size_t cases,
                  std::size_t partition)
      : m_model(model),
        m_acc(acc),
        m_threshold(threshold),
        m_nb_controls(controls),
        m_nb_cases(cases),
        m_part(partition)
    {
      m_counts.resize(m_nb_controls + m_nb_cases, 0);
    }

  public:
    void process(km::Kmer<KSIZE>& kmer, std::vector<count_type>& counts) override
    {
      Range<count_type> range_controls(counts, 0, m_nb_controls);
      Range<count_type> range_cases(counts, m_nb_controls, m_nb_cases);

      auto [p_value, sign, mean_ctr, mean_case] = m_model->process(range_controls, range_cases);

      spdlog::debug("P{}: {} {} {} {}", m_part, p_value, significance_to_char(sign), mean_ctr, mean_case);
      m_total++;

      if (p_value <= m_threshold)
      {

        km::Kmer<KSIZE> kmer_ = kmer;

        #ifndef WITH_POPSTRAT
          KmerSign<KSIZE> ks(std::move(kmer_), p_value, sign, mean_ctr, mean_case);
        #else
          for (std::size_t i=0; i < counts.size(); i++) { m_counts[i] = counts[i]; }
          KmerSign<KSIZE> ks(std::move(kmer_), p_value, sign, m_counts, mean_ctr, mean_case);
        #endif

        if (sign == Significance::CONTROL)
          m_sign_controls++;
        else
          m_sign_cases++;

        m_acc->push(std::move(ks));
        m_sign_kmer_per_part++;
      }
    }

    std::size_t total() const { return m_total; }
    std::size_t nb_sign() const { return m_sign_kmer_per_part; }

    std::tuple<std::size_t, std::size_t> nb_signs() const
    {
      return std::make_tuple(m_sign_controls, m_sign_cases);
    }

  protected:
    const std::shared_ptr<IModel<CMAX>> m_model {nullptr};
    std::size_t m_sign_kmer_per_part {0};
    std::size_t m_total {0};
    acc_t<KmerSign<KSIZE>> m_acc {nullptr};
    std::size_t m_nb_controls {0};
    std::size_t m_nb_cases {0};
    std::size_t m_part {0};
    double m_threshold {0};
    std::vector<double> m_counts;
    std::size_t m_sign_controls {0};
    std::size_t m_sign_cases {0};
};

template<std::size_t KSIZE, std::size_t CMAX>
class diff_observer_strat : public diff_observer<KSIZE, CMAX>
{
  using count_type = typename km::selectC<CMAX>::type;
  public:
    diff_observer_strat(
      const std::shared_ptr<IModel<CMAX>>& model,
      acc_t<KmerSign<KSIZE>> acc,
      double threshold,
      std::size_t controls,
      std::size_t cases,
      std::shared_ptr<Sampler<CMAX>> sampler,
      std::size_t partition
    ) : m_sampler(sampler), diff_observer<KSIZE, CMAX>(model, acc, threshold, controls, cases, partition) {}

    void process(km::Kmer<KSIZE>& kmer, std::vector<count_type>& counts) override
    {
      Range<count_type> range_controls(counts, 0, this->m_nb_controls);
      Range<count_type> range_cases(counts, this->m_nb_controls, this->m_nb_cases);

      m_sampler->sample(range_controls, range_cases);

      auto [p_value, sign, mean_ctr, mean_case] = this->m_model->process(range_controls, range_cases);
      spdlog::debug("P{}: {} {} {} {}", this->m_part, p_value, significance_to_char(sign), mean_ctr, mean_case);

      this->m_total++;

      if (p_value <= this->m_threshold)
      {
        km::Kmer<KSIZE> kmer_ = kmer;

        #ifndef WITH_POPSTRAT
          KmerSign<KSIZE> ks(std::move(kmer_), p_value, sign, mean_ctr, mean_case);
        #else
          for (std::size_t i=0; i < counts.size(); i++) { this->m_counts[i] = counts[i]; }
          KmerSign<KSIZE> ks(std::move(kmer_), p_value, sign, this->m_counts, mean_ctr, mean_case);
        #endif

        if (sign == Significance::CONTROL)
          this->m_sign_controls++;
        else
          this->m_sign_cases++;

        this->m_acc->push(std::move(ks));
        this->m_sign_kmer_per_part++;
      }
    }

  private:
    std::shared_ptr<Sampler<CMAX>> m_sampler {nullptr};
};

template<std::size_t KSIZE, std::size_t CMAX>
class global_merge
{
  public:
    global_merge(std::vector<std::vector<std::string>>& partition_paths,
                 std::vector<std::uint32_t>& ab_thresholds,
                 const std::shared_ptr<IModel<CMAX>>& model,
                 const std::vector<acc_t<KmerSign<KSIZE>>>& accumulators,
                 std::size_t kmer_size,
                 std::size_t nb_controls,
                 std::size_t nb_cases,
                 double threshold,
                 std::size_t nb_threads,
                 std::shared_ptr<Sampler<CMAX>> sampler)
      : m_part_paths(partition_paths),
        m_ab_thresholds(ab_thresholds),
        m_model(model),
        m_accs(accumulators),
        m_kmer_size(kmer_size),
        m_controls(nb_controls),
        m_cases(nb_cases),
        m_threshold(threshold),
        m_nb_threads(nb_threads),
        m_sampler(sampler)
    {}

    std::size_t merge()
    {
      ThreadPool pool(m_nb_threads);
      const std::size_t size = m_part_paths.size();

      std::vector<size_t> total_kmers(size);

      m_nb_signs.resize(size, 0);
      m_sign_controls.resize(size, 0);
      m_sign_cases.resize(size, 0);

      std::exception_ptr ep = nullptr;

      indicators::ProgressBar* pb = nullptr;

      if ((spdlog::get_level() != spdlog::level::debug) && isatty_stderr())
      {
        pb = get_progress_bar("progress", size, 50, indicators::Color::white, true);
        pb->print_progress();
      }

      for (std::size_t p = 0; p < size; p++)
      {
        auto partition_merger = [&ep, &total_kmers, p, pb, this](int id) {
          spdlog::debug("Process partition {}.", p);
          Timer mp_timer;

          km::KmerMerger<KSIZE, CMAX> km_merge(
            this->m_part_paths[p], this->m_ab_thresholds, this->m_kmer_size, 1, 0);

          km::imo_t<KSIZE, CMAX> diff {nullptr};

          if (!m_sampler)
            diff = std::make_shared<diff_observer<KSIZE, CMAX>>(
              this->m_model, this->m_accs[p], this->m_threshold,
              this->m_controls, this->m_cases, p);
          else
            diff = std::make_shared<diff_observer_strat<KSIZE, CMAX>>(
              this->m_model, this->m_accs[p], this->m_threshold,
              this->m_controls, this->m_cases, this->m_sampler, p);

          try { km_merge.merge(diff); } catch (...) { ep = std::current_exception(); }

          total_kmers[p] = dynamic_cast<diff_observer<KSIZE, CMAX>*>(diff.get())->total();
          this->m_nb_signs[p] = dynamic_cast<diff_observer<KSIZE, CMAX>*>(diff.get())->nb_sign();

          auto [co, ca] = dynamic_cast<diff_observer<KSIZE, CMAX>*>(diff.get())->nb_signs();

          this->m_sign_controls[p] += co;
          this->m_sign_cases[p] += ca;

          this->m_accs[p]->finish();

          spdlog::debug("Partition {} processed. ({})", p, mp_timer.formatted());
          if (pb)
            pb->tick();
        };

        pool.add_task(partition_merger);
      }

      pool.join_all();

      delete pb;

      if (ep != nullptr)
        rethrow_exception(ep);

      return std::accumulate(total_kmers.begin(), total_kmers.end(), 0ULL);
    }

    size_t nb_sign() const
    {
      return std::accumulate(m_nb_signs.begin(), m_nb_signs.end(), 0ULL);
    }

    std::tuple<size_t, size_t> signs() const
    {
      return std::make_tuple(
        std::accumulate(m_sign_controls.begin(), m_sign_controls.end(), 0ULL),
        std::accumulate(m_sign_cases.begin(), m_sign_cases.end(), 0ULL)
      );
    }

    private:
      std::vector<std::vector<std::string>>& m_part_paths;
      std::vector<std::uint32_t> m_ab_thresholds;
      const std::shared_ptr<IModel<CMAX>>& m_model;
      const std::vector<acc_t<KmerSign<KSIZE>>>& m_accs;
      std::size_t m_kmer_size {0};
      std::size_t m_controls {0};
      std::size_t m_cases {0};
      double m_threshold {0};
      std::size_t m_nb_threads {1};

      std::vector<size_t> m_nb_signs;
      std::vector<size_t> m_sign_controls;
      std::vector<size_t> m_sign_cases;

      std::shared_ptr<Sampler<CMAX>> m_sampler {nullptr};
};

};  // end of namespace kmdiff
