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
    diff_observer(const std::shared_ptr<Model<CMAX>>& model,
                  acc_t<KmerSign<KSIZE>> acc,
                  double threshold,
                  std::size_t controls,
                  std::size_t cases)
      : m_model(model),
        m_acc(acc),
        m_threshold(threshold),
        m_nb_controls(controls),
        m_nb_cases(cases)
    {}

  public:
    void process(km::Kmer<KSIZE>& kmer, std::vector<count_type>& counts)
    {
      Range<count_type> range_controls(counts, 0, m_nb_controls);
      Range<count_type> range_cases(counts, m_nb_controls, m_nb_cases);

      auto [p_value, sign, mean_ctr, mean_case] = m_model->process(range_controls, range_cases);

      m_total++;

      if (p_value <= m_threshold)
      {
        km::Kmer<KSIZE> kmer_ = kmer;
        KmerSign<KSIZE> ks(std::move(kmer_), p_value, sign, mean_ctr, mean_case);
        m_acc->push(std::move(ks));
        m_sign_kmer_per_part++;
      }
    }

    std::size_t total() const { return m_total; }
    std::size_t nb_sign() const { return m_sign_kmer_per_part; }

  private:
    const std::shared_ptr<Model<CMAX>> m_model {nullptr};
    std::size_t m_sign_kmer_per_part {0};
    std::size_t m_total {0};
    acc_t<KmerSign<KSIZE>> m_acc {nullptr};
    std::size_t m_nb_controls {0};
    std::size_t m_nb_cases {0};
    double m_threshold {0};
};

template<std::size_t KSIZE, std::size_t CMAX>
class global_merge
{
  public:
    global_merge(std::vector<std::vector<std::string>>& partition_paths,
                 std::vector<std::uint32_t>& ab_thresholds,
                 const std::shared_ptr<Model<CMAX>>& model,
                 const std::vector<acc_t<KmerSign<KSIZE>>>& accumulators,
                 std::size_t kmer_size,
                 std::size_t nb_controls,
                 std::size_t nb_cases,
                 double threshold,
                 std::size_t nb_threads)
      : m_part_paths(partition_paths),
        m_ab_thresholds(ab_thresholds),
        m_model(model),
        m_accs(accumulators),
        m_kmer_size(kmer_size),
        m_controls(nb_controls),
        m_cases(nb_cases),
        m_threshold(threshold),
        m_nb_threads(nb_threads)
    {}

    std::size_t merge()
    {
      ThreadPool pool(m_nb_threads);
      const std::size_t size = m_part_paths.size();

      std::vector<size_t> total_kmers(size);
      m_nb_signs.resize(size, 0);

      std::exception_ptr ep = nullptr;

      for (std::size_t p = 0; p < size; p++)
      {
        auto partition_merger = [&ep, &total_kmers, p, this](int id) {
          Timer mp_timer;

          spdlog::debug("Process partition {}...", p);

          km::KmerMerger<KSIZE, CMAX> km_merge(
            this->m_part_paths[p], this->m_ab_thresholds, this->m_kmer_size, 1, 0);

          km::imo_t<KSIZE, CMAX> diff = std::make_shared<diff_observer<KSIZE, CMAX>>(
            this->m_model, this->m_accs[p], this->m_threshold, this->m_controls, this->m_cases);

          try { km_merge.merge(diff); } catch (...) { ep = std::current_exception(); }

          total_kmers[p] = dynamic_cast<diff_observer<KSIZE, CMAX>*>(diff.get())->total();
          this->m_nb_signs[p] = dynamic_cast<diff_observer<KSIZE, CMAX>*>(diff.get())->nb_sign();
          this->m_accs[p]->finish();

          spdlog::debug("Partition {} processed", p);
        };

        pool.add_task(partition_merger);
      }
      pool.join_all();

      if (ep != nullptr)
        rethrow_exception(ep);

      return std::accumulate(total_kmers.begin(), total_kmers.end(), 0ULL);
    }

    size_t nb_sign() const
    {
      return std::accumulate(m_nb_signs.begin(), m_nb_signs.end(), 0ULL);
    }

    private:
      std::vector<std::vector<std::string>>& m_part_paths;
      std::vector<std::uint32_t> m_ab_thresholds;
      const std::shared_ptr<Model<CMAX>>& m_model;
      const std::vector<acc_t<KmerSign<KSIZE>>>& m_accs;
      std::size_t m_kmer_size;
      std::size_t m_controls;
      std::size_t m_cases;
      double m_threshold;
      std::size_t m_nb_threads;

      std::vector<size_t> m_nb_signs;
};

};  // end of namespace kmdiff
