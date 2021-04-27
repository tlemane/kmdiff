#include <kmdiff/validator.hpp>

namespace kmdiff {

Validator::Validator(const std::string& seq_path,
                     const std::string& kmer_path,
                     const std::string& out_path,
                     size_t kmer_size)
  : m_seq_path(seq_path), m_kmer_path(kmer_path), m_out_path(out_path), m_kmer_size(kmer_size)
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
  {
    if (record->core.tid >= 0)
    {
      targets[record->core.tid]++;
    }
  }
  nb_covered = std::count_if(targets.cbegin(), targets.cend(), [](int32_t v){
    return v > 0;
  });

  std::string stat_path = fmt::format("{}.infos", m_out_path);
  std::ofstream info(stat_path, std::ios::out);
  double sum = 0;
  for (size_t i = 0; i<nb_targets; i++)
  {
    size_t nbk = (header->target_len[i] - m_kmer_size + 1);
    double r = static_cast<double>(targets[i]) / static_cast<double>(nbk);
    info << header->target_name[i] << " " << std::to_string(nbk) << " ";
    info << std::to_string(targets[i]) << std::to_string(r) << "\n";
    sum += r;
  }
  spdlog::debug(sum / static_cast<double>(nb_targets));
}

}; // end of namespace kmdiff