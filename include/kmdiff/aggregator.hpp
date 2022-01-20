#pragma once

#include <string>
#include <vector>
#include <queue>
#include <mutex>

#include <spdlog/spdlog.h>
#include <kseq++/kseq++.hpp>
#include <kseq++/seqio.hpp>

#include <kmdiff/cmd/diff_opt.hpp>
#include <kmdiff/accumulator.hpp>
#include <kmdiff/kmer.hpp>
#include <kmdiff/kmtricks_utils.hpp>
#include <kmdiff/threadpool.hpp>
#include <kmdiff/blocking_queue.hpp>
#include <kmdiff/popstrat.hpp>
#include <kmdiff/kff_utils.hpp>
#include <kmdiff/progress.hpp>
#include <kmdiff/icorrector.hpp>

namespace kmdiff {
  using pb_t = indicators::ProgressBar*;

  template<size_t MAX_K>
  static void writer(BlockingQueue<KmerSign<MAX_K>>& queue,
                      const std::string& out_path,
                      const std::string& name,
                      bool kff,
                      kmtricks_config_t config,
                      std::size_t& count)
  {
    using seq_out_t = std::unique_ptr<klibpp::SeqStreamOut>;

    KmerSign<MAX_K> k;
    klibpp::KSeq record;
    seq_out_t out = nullptr;
    kff_w_t out_kff = nullptr;

    if (kff)
    {
      out_kff = std::make_unique<KffWriter>(out_path, config.kmer_size);
    }
    else
    {
      out = std::make_unique<klibpp::SeqStreamOut>(out_path.c_str());
      *out << klibpp::format::fasta;
    }

    while (queue.pop(k))
    {
      if (!kff)
      {
        record.name = fmt::format("{}_pval={:g}_control={}_case={}",
                                  count,
                                  k.m_pvalue,
                                  static_cast<std::size_t>(k.m_mean_control),
                                  k.m_mean_case);

        record.seq = k.m_kmer.to_string();
        *out << record;
      }
      else
      {
        out_kff->write(k);
      }
      count++;
    }
    if (out_kff) out_kff->close();
  }

  template<std::size_t KSIZE>
  class IAggregator
  {
    using ks_type = KmerSign<KSIZE>;

    public:
      IAggregator(std::vector<acc_t<ks_type>>& accumulators,
                  corrector_t corrector,
                  kmtricks_config_t config,
                  const std::string& output_dir,
                  bool kff,
                  std::size_t nb_threads,
                  pb_t pb = nullptr)
        : m_accumulators(accumulators),
          m_corrector(corrector),
          m_config(config),
          m_output(output_dir),
          m_kff(kff),
          m_nb_threads(nb_threads),
          m_pb(pb)
      {}

      virtual void run() = 0;

      std::tuple<std::size_t, std::size_t> counts() const
      {
        return std::make_tuple(m_control_count, m_case_count);
      }

    protected:
      std::vector<acc_t<ks_type>>& m_accumulators;
      corrector_t m_corrector {nullptr};
      kmtricks_config_t m_config;

      const std::string& m_output;

      std::size_t m_control_count {0};
      std::size_t m_case_count {0};
      std::size_t m_nb_significant {0};

      bool m_kff {false};
      std::size_t m_nb_threads {1};

      pb_t m_pb {nullptr};
  };

  template<std::size_t KSIZE>
  class aggregator : public IAggregator<KSIZE>
  {
    using ks_type = KmerSign<KSIZE>;

    public:
      aggregator(std::vector<acc_t<ks_type>>& accumulators,
                 corrector_t corrector,
                 kmtricks_config_t config,
                 const std::string& output_dir,
                 bool kff,
                 std::size_t nb_threads,
                 pb_t pb = nullptr)
        : IAggregator<KSIZE>(accumulators, corrector, config, output_dir, kff, nb_threads, pb)
      {}

      static void worker(BlockingQueue<KmerSign<KSIZE>>& controls_queue,
                         BlockingQueue<KmerSign<KSIZE>>& cases_queue,
                         acc_t<ks_type>& accumulator,
                         corrector_t corrector,
                         kmtricks_config_t config,
                         std::size_t partition,
                         pb_t pb,
                         int thread_id)
      {
        while (auto& o = accumulator->get())
        {
          auto& kref = o.value();

          bool keep = true;
          keep = corrector->apply(kref.m_pvalue);

          if (keep)
          {
            if (kref.m_sign == Significance::CONTROL)
            {
              controls_queue.push(std::move(kref));
            }
            else
            {
              cases_queue.push(std::move(kref));
            }
          }
        }

        controls_queue.end_signal(partition);
        cases_queue.end_signal(partition);

        if (pb)
          pb->tick();
      }

      void run() override
      {
        const auto& nb_part = this->m_config.nb_partitions;
        const auto& nb_threads = this->m_nb_threads;

        BlockingQueue<ks_type> cases_queue(50000, nb_part);
        BlockingQueue<ks_type> controls_queue(50000, nb_part);

        ThreadPool pool(nb_threads < 2 ? 1 : nb_threads);

        for (std::size_t p = 0; p < nb_part; p++)
        {
          auto task = std::bind(&aggregator<KSIZE>::worker,
                                std::ref(controls_queue),
                                std::ref(cases_queue),
                                std::ref(this->m_accumulators[p]),
                                this->m_corrector,
                                this->m_config,
                                p,
                                this->m_pb,
                                std::placeholders::_1);
          pool.add_task(task);
        }

        std::string ext = this->m_kff ? ".kff" : ".fasta";
        std::string control_out = fmt::format("{}/control_kmers{}", this->m_output, ext);
        std::string case_out = fmt::format("{}/case_kmers{}", this->m_output, ext);

        auto control_writer = std::thread(&writer<KSIZE>,
                                          std::ref(controls_queue),
                                          control_out,
                                          "control",
                                          this->m_kff,
                                          this->m_config,
                                          std::ref(this->m_control_count));

        auto case_writer = std::thread(&writer<KSIZE>,
                                       std::ref(cases_queue),
                                       case_out,
                                       "case",
                                       this->m_kff,
                                       this->m_config,
                                       std::ref(this->m_case_count));

        control_writer.join();
        case_writer.join();
        pool.join_all();
      }
  };

  template<std::size_t KSIZE>
  class sorted_aggregator : public IAggregator<KSIZE>
  {
    using ks_type = KmerSign<KSIZE>;

    public:
      sorted_aggregator(std::vector<acc_t<ks_type>>& accumulators,
                        corrector_t corrector,
                        kmtricks_config_t config,
                        const std::string& output_dir,
                        bool kff,
                        std::size_t nb_threads,
                        pb_t pb = nullptr)
        : IAggregator<KSIZE>(accumulators, corrector, config, output_dir, kff, nb_threads, pb)
      {}


      void run() override
      {
        const auto& nb_part = this->m_config.nb_partitions;
        const auto& nb_threads = this->m_nb_threads;

        BlockingQueue<ks_type> cases_queue(50000, nb_part);
        BlockingQueue<ks_type> controls_queue(50000, nb_part);

        ThreadPool pool(nb_threads < 2 ? 1 : nb_threads);

        for (std::size_t p = 0; p < nb_part; p++)
        {
          auto task = std::bind(&sorted_aggregator<KSIZE>::worker,
                                std::ref(m_pqueue),
                                std::ref(this->m_accumulators[p]),
                                p,
                                std::ref(m_lock),
                                std::placeholders::_1);
          pool.add_task(task);
        }
        pool.join_all();

        std::string ext = this->m_kff ? ".kff" : ".fasta";
        std::string control_out = fmt::format("{}/control_kmers{}", this->m_output, ext);
        std::string case_out = fmt::format("{}/case_kmers{}", this->m_output, ext);

        auto control_writer = std::thread(&writer<KSIZE>,
                                          std::ref(controls_queue),
                                          control_out,
                                          "control",
                                          this->m_kff,
                                          this->m_config,
                                          std::ref(this->m_control_count));

        auto case_writer = std::thread(&writer<KSIZE>,
                                       std::ref(cases_queue),
                                       case_out,
                                       "case",
                                       this->m_kff,
                                       this->m_config,
                                       std::ref(this->m_case_count));

        std::size_t c = m_pqueue.size() / nb_part;

        std::size_t i = 0;

        while (!m_pqueue.empty())
        {
          auto ks = m_pqueue.top();

          if (!this->m_corrector->apply(ks.m_pvalue))
            break;

          if (ks.m_sign == Significance::CONTROL)
          {
            controls_queue.push(std::move(ks));
          }
          else
          {
            cases_queue.push(std::move(ks));
          }

          if (this->m_pb)
          {
            if ((++i % c) == 0)
            {
              this->m_pb->tick();
            }
          }
          m_pqueue.pop();
        }

        if (this->m_pb && !this->m_pb->is_completed())
        {
          this->m_pb->set_progress(nb_part);
        }

        controls_queue.end_signal();
        cases_queue.end_signal();

        control_writer.join();
        case_writer.join();
      }

    private:
      static void worker(std::priority_queue<KmerSign<KSIZE>>& pq,
                         acc_t<ks_type>& accumulator,
                         std::size_t partition,
                         spinlock& lock,
                         int thread_id)
      {
        while (auto& o = accumulator->get())
        {
          std::unique_lock<spinlock> lock_(lock);
          pq.push(o.value());
        }
      }

    private:
      std::priority_queue<KmerSign<KSIZE>> m_pqueue;
      spinlock m_lock{};
  };

}; // end of namespace kmdiff

