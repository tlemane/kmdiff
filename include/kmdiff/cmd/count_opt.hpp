#pragma once

#include <kmdiff/cmd/cmd_common.hpp>

namespace kmdiff {

  struct count_options : kmdiff_options
  {
    std::string file;
    std::string dir;
    int kmer_size;
    int abundance_min;
    int recurrence_min;
    int memory;
    int minimizer_type;
    int minimizer_size;
    int repartition_type;
    int nb_partitions;

    std::string display()
    {
      std::stringstream ss;
      ss << this->global_display();
      KRECORD(ss, file);
      KRECORD(ss, dir);
      KRECORD(ss, kmer_size);
      KRECORD(ss, abundance_min);
      KRECORD(ss, recurrence_min);
      KRECORD(ss, memory);
      KRECORD(ss, minimizer_type);
      KRECORD(ss, minimizer_size);
      KRECORD(ss, repartition_type);
      KRECORD(ss, nb_partitions);
      return ss.str();
    }
  };

  using count_options_t = std::shared_ptr<struct count_options>;

} // end of namespace kmdiff
