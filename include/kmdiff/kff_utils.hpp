
// std
#include <cstdint>
#include <string>
#include <memory>

// ext
#include <kff_io.hpp>

// int
#include <kmdiff/kmer.hpp>

namespace kmdiff {

using kff_t = std::shared_ptr<Kff_file>;
using kff_raw_t = std::shared_ptr<Section_Raw>;

class KffWriter
{
public:
  KffWriter(const std::string& path, size_t kmer_size)
    : m_kmer_size(kmer_size)
  {
    m_kff_file = std::make_shared<Kff_file>(path, "w");
    Section_GV sgv(m_kff_file.get());
    sgv.write_var("k", m_kmer_size);
    sgv.write_var("max", 1);
    sgv.write_var("data_size", 0);
    sgv.close();
    m_kff_sec = std::make_shared<Section_Raw>(m_kff_file.get());
  }

  template<size_t MAX_K>
  void write(KmerSign<MAX_K>& kmer)
  {
    uint8_t encoded[1024];
    m_kmer = kmer.m_kmer.to_string();
    encode_sequence(encoded);
    m_kff_sec->write_compacted_sequence(encoded, m_kmer_size, nullptr);
  }

  void close()
  {
    m_kff_sec->close();
    m_kff_file->close();
  }

  uint8_t uint8_packing(std::string sequence)
  {
    size_t size = sequence.size();
    uint8_t val = 0;
    for (size_t i = 0; i < size; i++)
    {
      val <<= 2;
      val += (sequence[i] >> 1) & 0b11;
    }
    return val;
  }
  
  void encode_sequence(uint8_t* encoded)
  {
    size_t size = m_kmer.length();
    size_t remnant = size % 4;
    if (remnant > 0)
    {
      encoded[0] = uint8_packing(m_kmer.substr(0, remnant));
      encoded += 1;
    }

    size_t nb_uint_needed = size / 4;
    for (size_t i = 0; i < nb_uint_needed; i++)
    {
      encoded[i] = uint8_packing(m_kmer.substr(remnant + 4 * i, 4));
    }
  }
private:
  kff_t m_kff_file {nullptr};
  kff_raw_t m_kff_sec {nullptr};
  size_t m_kmer_size;
  std::string m_kmer;
};

using kff_w_t = std::unique_ptr<KffWriter>;

}; // end of namespace kmdiff