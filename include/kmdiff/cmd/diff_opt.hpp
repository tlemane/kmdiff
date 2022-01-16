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
    CorrectionType correction;
    bool in_memory;
    bool kff;

  #ifdef WITH_POPSTRAT
    bool pop_correction;
    double kmer_pca;
    size_t ploidy;
    bool is_diploid;
    size_t npc;
    std::string covariates;
  #endif

  #ifdef KMDIFF_DEV_MODE
    double learning_rate;
    size_t max_iteration;
    double epsilon;
    bool stand {false};
  #endif

    std::string display()
    {
      std::stringstream ss;
      ss << this->global_display();
      KRECORD(ss, kmtricks_dir);
      KRECORD(ss, output_directory);
      KRECORD(ss, nb_controls);
      KRECORD(ss, nb_cases);
      KRECORD(ss, threshold);
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
  #endif
  #ifdef KMDIFF_DEV_MODE
      KRECORD(ss, learning_rate);
      KRECORD(ss, max_iteration);
      KRECORD(ss, epsilon);
      KRECORD(ss, stand)
  #endif
      return ss.str();
    }
  };

using diff_options_t = std::shared_ptr<struct diff_options>;

}
