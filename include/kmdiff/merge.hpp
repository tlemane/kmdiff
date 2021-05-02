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

// ext
#define _KM_LIB_INCLUDE_
#include <kmtricks/io.hpp>
#include <kmtricks/merger.hpp>
#include <kmtricks/utilities.hpp>

// int
#include <kmdiff/accumulator.hpp>
#include <kmdiff/model.hpp>
#include <kmdiff/threadpool.hpp>
#include <kmdiff/utils.hpp>
#include <kmdiff/blocking_queue.hpp>
#include <kmdiff/popstrat.hpp>

namespace kmdiff
{
template <size_t MAX_KMER_SIZE, size_t MAX_COUNT>
class GlobalMerge
{
  using ktype = typename selectK<MAX_KMER_SIZE>::type;
  using ctype = typename selectC<MAX_COUNT>::type;

 public:
  GlobalMerge() = delete;
  GlobalMerge(const GlobalMerge&) = delete;
  GlobalMerge& operator=(const GlobalMerge&) = delete;
  GlobalMerge(const GlobalMerge&&) = delete;
  GlobalMerge& operator=(const GlobalMerge&&) = delete;

  GlobalMerge(
      std::vector<std::string>& fofs,
      std::vector<uint32_t>& a_min,
      const int kmer_size,
      const size_t nb_controls,
      const size_t nb_cases,
      double threshold,
      const std::string& output_directory,
      const int threads,
      const std::shared_ptr<Model<MAX_COUNT>>& model,
      const std::vector<acc_t<KmerSign<MAX_KMER_SIZE>>>& accumulators,
      const std::string& path_geno,
      const std::string& path_snp,
      bool with_popstrat,
      double for_pca)
      : m_fofs(fofs),
        m_a_min(a_min),
        m_kmer_size(kmer_size),
        m_output_directory(output_directory),
        m_threads(threads),
        m_nb_controls(nb_controls),
        m_nb_cases(nb_cases),
        m_model(model),
        m_threshold(threshold),
        m_accumulators(accumulators),
        m_path_geno(path_geno),
        m_path_snp(path_snp),
        m_with_popstrat(with_popstrat),
        m_for_pca(for_pca)
  {
    m_total_kmers.resize(m_fofs.size(), 0);
    m_sign_per_part.resize(m_fofs.size(), 0);
  }

  size_t merge()
  {
    ThreadPool pool(m_threads);
    const size_t size = m_fofs.size();

    std::exception_ptr ep = nullptr;

    std::shared_ptr<EigGenoFile<MAX_COUNT>> geno_file = nullptr;
    std::shared_ptr<EigSnpFile> snp_file = nullptr;

#ifdef WITH_POPSTRAT
    if (m_with_popstrat)
    {
      geno_file = std::make_shared<EigGenoFile<MAX_COUNT>>(m_path_geno);
      snp_file = std::make_shared<EigSnpFile>(m_path_snp);
    }
#endif

    for (int f = 0; f < m_fofs.size(); f++)
    {
      auto partition_merger = [&size, &ep, this, f, geno_file, snp_file](int id) {
        Timer partition_timer;
        spdlog::debug("Process partition {}...", f);

        try
        {
          km::Merger<ktype, ctype, km::KmerFile<km::IN, ktype, ctype>> merger(
              this->m_fofs[f], this->m_a_min, 1, 0, false, 1, false);

          Range<ctype> r_controls(merger.counts, 0, m_nb_controls);
          Range<ctype> r_cases(merger.counts, m_nb_controls, m_nb_cases);

          while (!merger.end)
          {
            merger.next();
            if (merger.keep)
            {
#ifdef WITH_POPSTRAT
              if (this->m_with_popstrat && (m_dist(m_generator) < m_for_pca))
              {
                geno_file->push(r_controls, r_cases);
                snp_file->push();
              }
#endif
              this->m_total_kmers[f]++;

              auto [p_value, sign, mean_ctr, mean_case] = m_model->process(r_controls, r_cases);

              Kmer<MAX_KMER_SIZE> k;
              k.from_km(&merger.m_khash, m_kmer_size);
              if (p_value <= (this->m_threshold / static_cast<double>(100000)))
              {
                KmerSign<MAX_KMER_SIZE> ks(std::move(k), p_value, sign, mean_ctr, mean_case);
                this->m_accumulators[f]->push(std::move(ks));
                this->m_sign_per_part[f]++;
              }
            }
          }
          this->m_accumulators[f]->finish();
          spdlog::debug(
              "Partition {} processed. ({} seconds)", f,
              partition_timer.elapsed<std::chrono::seconds>().count());
        }
        catch (...)
        {
          ep = std::current_exception();
        }
      };
      pool.add_task(partition_merger);
    }
    pool.join_all();

    if (ep != nullptr) rethrow_exception(ep);

    return std::accumulate(m_total_kmers.begin(), m_total_kmers.end(), 0ULL);
  }

  size_t nb_sign()
  {
    return std::accumulate(m_sign_per_part.begin(), m_sign_per_part.end(), 0ULL);
  }

 private:
  std::vector<std::string>& m_fofs;
  std::vector<uint32_t>& m_a_min;
  std::string m_output_directory;
  std::vector<size_t> m_total_kmers;
  std::vector<size_t> m_sign_per_part;
  size_t m_nb_controls{0};
  size_t m_nb_cases{0};
  int m_threads{0};
  int m_kmer_size{0};
  double m_threshold{0};
  const std::shared_ptr<Model<MAX_COUNT>> m_model;
  const std::vector<acc_t<KmerSign<MAX_KMER_SIZE>>>& m_accumulators;
  std::string m_path_geno;
  std::string m_path_snp;
  bool m_with_popstrat;
  std::default_random_engine m_generator;
  std::uniform_real_distribution<double> m_dist {0.0, 1.0};
  double m_for_pca;
};

};  // end of namespace kmdiff
