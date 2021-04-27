#include <kmdiff/validator.hpp>

namespace kmdiff {

Validator::Validator(const std::string& seq_path,
                     const std::string& kmer_path,
                     const std::string& out_path)
  : m_seq_path(seq_path), m_kmer_path(kmer_path), m_out_path(out_path)
{
}

void Validator::align(size_t seed_size, size_t nb_threads)
{
  std::string bbmap_bin = command_exists(get_binary_dir(), "bbmap.sh");
  std::string bbmap_align = "ref={} in={} out={} k={} threads={} nodisk";
  std::string cmd = fmt::format(bbmap_align,
                                m_seq_path,
                                m_kmer_path,
                                m_out_path,
                                seed_size,
                                nb_threads);
  exec_external_cmd(bbmap_bin, cmd);
}

void Validator::valid(size_t& nb_targets, size_t& nb_covered)
{
  BAMReader reader(m_out_path);
  bam_hdr_t* header = reader.get_header();

  nb_targets = header->n_targets;
  bam1_t* record;

  std::vector<uint8_t> targets(nb_targets, 0);
  while (record = reader.next())
    if (record->core.tid >= 0)
      targets[record->core.tid] = 1;

  nb_covered = std::accumulate(targets.begin(), targets.end(), 0ULL);
}

}; // end of namespace kmdiff