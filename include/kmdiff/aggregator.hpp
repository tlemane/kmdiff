#pragma once
// std
#include <string>
#include <vector>
#include <queue>
#include <mutex>

// ext
#include <spdlog/spdlog.h>
#include <kseq++/kseq++.hpp>
#include <kseq++/seqio.hpp>

// int
#include <kmdiff/cmd/diff.hpp>
#include <kmdiff/merge.hpp>
#include <kmdiff/model.hpp>
#include <kmdiff/accumulator.hpp>
#include <kmdiff/kmer.hpp>
#include <kmdiff/kmtricks_utils.hpp>
#include <kmdiff/cmd/diff.hpp>
#include <kmdiff/threadpool.hpp>
#include <kmdiff/blocking_queue.hpp>


namespace kmdiff {

template<size_t MAX_K>
class IAggregator
{
public:
  IAggregator(std::vector<acc_t<KmerSign<MAX_K>>>& accumulators,
              std::shared_ptr<ICorrector> corrector,
              diff_options_t opt,
              kmtricks_config_t config)
    : m_accumulators(accumulators), m_corrector(corrector), m_opt(opt), m_config(config)
  {
  }

  virtual void run() = 0;

protected:
  size_t m_nb_significant;
  std::vector<acc_t<KmerSign<MAX_K>>>& m_accumulators;
  std::shared_ptr<ICorrector> m_corrector;
  diff_options_t m_opt;
  kmtricks_config_t m_config;
};

template<size_t MAX_K>
static void writer(BlockingQueue<KmerSign<MAX_K>>& queue,
                    const std::string& out_path, const std::string& name,
                    diff_options_t opt, kmtricks_config_t config)
{
  using seq_out_t = std::unique_ptr<klibpp::SeqStreamOut>;

  KmerSign<MAX_K> k;
  klibpp::KSeq record;
  seq_out_t out = nullptr;
  kff_w_t out_kff = nullptr;

  if (opt->kff)
    out_kff = std::make_unique<KffWriter>(out_path, config.kmer_size);
  else
    out = std::make_unique<klibpp::SeqStreamOut>(out_path.c_str());

  size_t count;
  while (queue.pop(k))
  {
    if (!opt->kff)
    {
      record.name = fmt::format("{}_pval={}_control={}_case={}",
                                count, k.m_pvalue, k.m_mean_control, k.m_mean_case);
      record.seq = k.m_kmer.to_string();
      *out << klibpp::format::fasta << record;
    }

    else
    {
      out_kff->write(k);
    }
    count++;
  }
  if (out_kff) out_kff->close();
  spdlog::info("Over-represented k-mers ({}) in {} dumped at {}", count, name, out_path);
}

template<size_t MAX_K>
class BonferonniAggregator : public IAggregator<MAX_K>
{
public:
  BonferonniAggregator(std::vector<acc_t<KmerSign<MAX_K>>>& accumulators,
                       std::shared_ptr<ICorrector> corrector,
                       diff_options_t opt,
                       kmtricks_config_t config)
    : IAggregator<MAX_K>(accumulators, corrector, opt, config)
  {}

  void run()
  {
    BlockingQueue<KmerSign<MAX_K>> case_queue(50000, this->m_config.nb_partitions);
    BlockingQueue<KmerSign<MAX_K>> control_queue(50000, this->m_config.nb_partitions);

    ThreadPool pool(this->m_opt->nb_threads < 2 ? 1 : this->m_opt->nb_threads);
    for (size_t p = 0; p < this->m_config.nb_partitions; p++)
    {
      auto task = std::bind(&BonferonniAggregator<MAX_K>::worker,
                            std::ref(control_queue), std::ref(case_queue),
                            std::ref(this->m_accumulators), this->m_corrector,
                            this->m_opt, this->m_config, p, std::placeholders::_1);
      pool.add_task(task);
    }

    std::string ext = this->m_opt->kff ? ".kff" : ".fasta";
    std::string control_out = fmt::format("{}/control_kmers{}", this->m_opt->output_directory, ext);
    std::string case_out = fmt::format("{}/case_kmers{}", this->m_opt->output_directory, ext);

    auto control_writer = std::thread(&writer<MAX_K>,
                                      std::ref(control_queue), control_out, "control",
                                      this->m_opt, this->m_config);
    auto case_writer = std::thread(&writer<MAX_K>,
                                   std::ref(case_queue), case_out, "case",
                                   this->m_opt, this->m_config);
    control_writer.join();
    case_writer.join();
    pool.join_all();
  }

  static void worker(BlockingQueue<KmerSign<MAX_K>>& control_queue,
                     BlockingQueue<KmerSign<MAX_K>>& case_queue,
                     std::vector<acc_t<KmerSign<MAX_K>>>& accumulators,
                     std::shared_ptr<ICorrector> corrector,
                     diff_options_t opt, kmtricks_config_t config,
                     size_t partition, int thread_id)
  {
    while (std::optional<KmerSign<MAX_K>>& o = accumulators[partition]->get())
    {
      bool keep = true;
      if (corrector && (opt->correction == CorrectionType::BONFERRONI))
        keep = corrector->apply((*o).m_pvalue);
      if (keep)
      {
        if ((*o).m_sign == Significance::CONTROL)
          control_queue.push(std::move(*o));
        else if ((*o).m_sign == Significance::CASE)
          case_queue.push(std::move(*o));
      }
    }
    control_queue.end_signal(partition);
    case_queue.end_signal(partition);
  }
};

template<size_t MAX_K>
class BenjaminiAggregator : public IAggregator<MAX_K>
{
public:
  BenjaminiAggregator(std::vector<acc_t<KmerSign<MAX_K>>>& accumulators,
                      std::shared_ptr<ICorrector> corrector,
                      diff_options_t opt,
                      kmtricks_config_t config)
    : IAggregator<MAX_K>(accumulators, corrector, opt, config)
  {}

  std::priority_queue<KmerSign<MAX_K>> m_pqueue;

  void run() override
  {
    BlockingQueue<KmerSign<MAX_K>> case_queue(50000, this->m_config.nb_partitions);
    BlockingQueue<KmerSign<MAX_K>> control_queue(50000, this->m_config.nb_partitions);

    ThreadPool pool(this->m_opt->nb_threads < 2 ? 1 : this->m_opt->nb_threads);
    for (size_t p = 0; p < this->m_config.nb_partitions; p++)
    {
      auto task = std::bind(&BenjaminiAggregator<MAX_K>::worker,
                            std::ref(m_pqueue),
                            std::ref(this->m_accumulators),
                            p, std::ref(m_mutex),
                            std::placeholders::_1);
      pool.add_task(task);
    }
    pool.join_all();

    std::string ext = this->m_opt->kff ? ".kff" : ".fasta";
    std::string control_out = fmt::format("{}/control_kmers{}", this->m_opt->output_directory, ext);
    std::string case_out = fmt::format("{}/case_kmers{}", this->m_opt->output_directory, ext);

    auto control_writer = std::thread(&writer<MAX_K>,
                                      std::ref(control_queue), control_out, "control",
                                      this->m_opt, this->m_config);
    auto case_writer = std::thread(&writer<MAX_K>,
                                   std::ref(case_queue), case_out, "case",
                                   this->m_opt, this->m_config);

    while (!m_pqueue.empty())
    {
      KmerSign<MAX_K> k = std::move(m_pqueue.top());
      if (!this->m_corrector->apply(k.m_pvalue))
        break;
      if (k.m_sign == Significance::CONTROL)
        control_queue.push(std::move(k));
      else if (k.m_sign == Significance::CASE)
        case_queue.push(std::move(k));
      m_pqueue.pop();
    }

    control_queue.end_signal();
    case_queue.end_signal();

    control_writer.join();
    case_writer.join();
  }

private:
  static void worker(std::priority_queue<KmerSign<MAX_K>>& pq,
                     std::vector<acc_t<KmerSign<MAX_K>>>& accumulators,
                     size_t partition, std::mutex& mutex, int thread_id)
  {
    while (std::optional<KmerSign<MAX_K>>& o = accumulators[partition]->get())
    {
      std::unique_lock<std::mutex> lock(mutex);
      pq.push((*o));
    }
  }

private:
std::mutex m_mutex;
};

template<size_t MAX_K>
auto make_uncorrected_aggregator(std::vector<acc_t<KmerSign<MAX_K>>>& accumulators,
                                 diff_options_t opt,
                                 kmtricks_config_t config)
{
  return std::make_unique<BonferonniAggregator<MAX_K>>(accumulators,
                                                       nullptr,
                                                       opt,
                                                       config);
}


}; // end of namespace kmdiff