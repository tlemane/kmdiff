#include <kmdiff/reference.hpp>
#include <kmdiff/simulator.hpp>

namespace kmdiff {

Reference::Reference(const std::string& path)
  : m_path(path)
{
  int status = fai_build(m_path.c_str());
  m_faidx = fai_load_format(m_path.c_str(), FAI_FASTA);
}

Reference::~Reference()
{
  fai_destroy(m_faidx);
}

std::string Reference::fetch(const std::string& name, size_t start, size_t end)
{
  std::string reg = fmt::format("{}:{}-{}", name, start, end);
  int len;
  char* seq = fai_fetch(m_faidx, reg.c_str(), &len);
  
  std::string ret(seq);
  delete seq;
  return ret;
}

std::tuple<std::string, std::string> Reference::get_kmer_view(SV& sv, size_t kmer_size)
{
  std::string ref_seq;
  std::string mut_seq;
  
  if (sv.m_type == visor_map.at("DEL"))
  {
    ref_seq = fetch(sv.m_chr, sv.m_start - kmer_size + 1, sv.m_end + kmer_size - 1);
    mut_seq = fetch(sv.m_chr, sv.m_start - kmer_size + 1, sv.m_start - 1) +
              fetch(sv.m_chr, sv.m_end + 1, sv.m_end + kmer_size - 1);
  }
  else if (sv.m_type == visor_map.at("INS"))
  {
    ref_seq = fetch(sv.m_chr, sv.m_end - kmer_size + 2, sv.m_end + kmer_size - 1);
    mut_seq = fetch(sv.m_chr, sv.m_end - kmer_size + 2, sv.m_end) +
              sv.m_seq +
              fetch(sv.m_chr, sv.m_end + 1, sv.m_end + kmer_size - 1);
  }
  else if (sv.m_type == visor_map.at("INV"))
  {
    ref_seq = fetch(sv.m_chr, sv.m_start - kmer_size + 1, sv.m_end + kmer_size - 1);
    mut_seq = fetch(sv.m_chr, sv.m_start - kmer_size + 1, sv.m_start - 1) +
              str_rev_comp(fetch(sv.m_chr, sv.m_start, sv.m_end)) +
              fetch(sv.m_chr, sv.m_end + 1, sv.m_end + kmer_size - 1);
  }
  else if (sv.m_type == visor_map.at("DUP"))
  {
    int nb = std::stoi(sv.m_seq);
    std::string v = fetch(sv.m_chr, sv.m_start, sv.m_end);
    std::stringstream ss;
    for (int i=0; i<nb; i++) ss << v;
  
    ref_seq = fetch(sv.m_chr, sv.m_start - kmer_size + 1, sv.m_end + kmer_size - 1);
    mut_seq = fetch(sv.m_chr, sv.m_start - kmer_size + 1, sv.m_start - 1) +
              ss.str() +
              fetch(sv.m_chr, sv.m_end + 1, sv.m_end + kmer_size - 1);
  }
  else if (sv.m_type == visor_map.at("IDUP"))
  {
    int nb = std::stoi(sv.m_seq);
    std::string v = fetch(sv.m_chr, sv.m_start, sv.m_end);
    std::stringstream ss; ss << v;
    for (int i=1; i<nb; i++) ss << str_rev_comp(v);
  
    ref_seq = fetch(sv.m_chr, sv.m_start - kmer_size + 1, sv.m_end + kmer_size - 1);
    mut_seq = fetch(sv.m_chr, sv.m_start - kmer_size + 1, sv.m_start - 1) +
              ss.str() +
              fetch(sv.m_chr, sv.m_end + 1, sv.m_end + kmer_size - 1);
  }
  else
  {
    throw std::runtime_error(fmt::format("{} not supported.", sv.m_type));
  }
  return std::make_tuple(ref_seq, mut_seq);
}

};