#pragma once

#include <kmdiff/cmd/cmd_common.hpp>

namespace kmdiff {
  struct diff_options : kmdiff_options
  {
    std::string kmtricks_dir;
    std::string output_directory;
    size_t nb_controls;
    size_t nb_cases;
    double threshold;
    double cutoff;
    CorrectionType correction;
    bool in_memory;
    bool cpr;
    bool kff;

    std::string model_lib_path;
    std::string model_config;

    bool pop_correction;
    double kmer_pca;
    size_t ploidy;
    bool is_diploid;
    size_t npc;
    std::string covariates;
    std::string gender;

    double learning_rate;
    size_t max_iteration;
    double epsilon;
    bool stand {false};
    bool irls {false};

    bool keep_tmp {false};

    std::size_t seed;
    std::size_t log_size;

    std::size_t total_kmers;

    std::string display()
    {
      std::stringstream ss;
      ss << this->global_display();
      KRECORD(ss, kmtricks_dir);
      KRECORD(ss, output_directory);
      KRECORD(ss, nb_controls);
      KRECORD(ss, nb_cases);
      KRECORD(ss, threshold);
      KRECORD(ss, cutoff);
      KRECORD(ss, correction_type_str(correction));
      KRECORD(ss, in_memory);
      KRECORD(ss, kff);
  #ifdef WITH_POPSTRAT
      KRECORD(ss, pop_correction);
      KRECORD(ss, kmer_pca);
      KRECORD(ss, ploidy);
      KRECORD(ss, is_diploid);
      KRECORD(ss, npc);
      KRECORD(ss, covariates);
      KRECORD(ss, gender);
  #endif
      KRECORD(ss, learning_rate);
      KRECORD(ss, max_iteration);
      KRECORD(ss, epsilon);
      KRECORD(ss, stand);
      KRECORD(ss, seed);
      return ss.str();
    }
  };

using diff_options_t = std::shared_ptr<struct diff_options>;

inline void dump_opt(diff_options_t opt, const std::string& path)
{
  std::ofstream out(path, std::ios::out | std::ios::binary);

  out.write(reinterpret_cast<char*>(&opt->threshold), sizeof(opt->threshold));
  out.write(reinterpret_cast<char*>(&opt->cutoff), sizeof(opt->cutoff));
  out.write(reinterpret_cast<char*>(&opt->correction), sizeof(opt->correction));
  out.write(reinterpret_cast<char*>(&opt->pop_correction), sizeof(opt->pop_correction));
  out.write(reinterpret_cast<char*>(&opt->kmer_pca), sizeof(opt->kmer_pca));
  out.write(reinterpret_cast<char*>(&opt->npc), sizeof(opt->npc));
}

inline diff_options_t load_opt(const std::string& path)
{
  std::ifstream in(path, std::ios::out | std::ios::binary);

  auto opt = std::make_shared<struct diff_options>();

  in.read(reinterpret_cast<char*>(&opt->threshold), sizeof(opt->threshold));
  in.read(reinterpret_cast<char*>(&opt->cutoff), sizeof(opt->cutoff));
  in.read(reinterpret_cast<char*>(&opt->correction), sizeof(opt->correction));
  in.read(reinterpret_cast<char*>(&opt->pop_correction), sizeof(opt->pop_correction));
  in.read(reinterpret_cast<char*>(&opt->kmer_pca), sizeof(opt->kmer_pca));
  in.read(reinterpret_cast<char*>(&opt->npc), sizeof(opt->npc));

  return opt;
}

inline unsigned compare_opt(diff_options_t opt, diff_options_t prev)
{
  unsigned r = 0;

  if (opt->threshold != prev->threshold || opt->cutoff != prev->cutoff)
    r |= 0b1;

  if (prev->pop_correction && opt->pop_correction)
  {
    if (opt->kmer_pca != prev->kmer_pca)
      r |= 0b11;
    if (opt->npc != prev->npc)
      r |= 0b10;
  }

  if (!prev->pop_correction && opt->pop_correction)
  {
    r |= 0b11; // 0b1 for sampling
  }

  if (opt->correction != prev->correction)
    r |= 0b100;

  if (prev->pop_correction && !opt->pop_correction)
    r |= 0b100;

  return r;
}

}
