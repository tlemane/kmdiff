// std
#include <string>
#include <vector>
#include <tuple>

// ext
#include <robin_hood.h>

// int
#include <kmdiff/io/bam.hpp>
#include <kmdiff/utils.hpp>

namespace kmdiff {

class Validator
{
public:
  Validator(const std::string& seq_path,
            const std::string& kmer_path,
            const std::string& out_path,
            size_t kmer_size);

  void align(size_t seed_size, size_t nb_threads);
  void valid(size_t& nb_targets, size_t& nb_covered);

private:
  std::string m_seq_path;
  std::string m_kmer_path;
  std::string m_out_path;
  size_t m_kmer_size;
};

}; // end of namespace kmdiff