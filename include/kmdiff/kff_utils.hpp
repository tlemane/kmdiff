
// std
#include <cstdint>
#include <string>
#include <memory>

// ext
#include <kff_io.hpp>

// int
#include <kmdiff/kmer.hpp>

namespace kmdiff {

inline uint8_t uint8_packing(std::string sequence)
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

inline void encode_sequence(std::string sequence, uint8_t* encoded)
{
  size_t size = sequence.length();
  size_t remnant = size % 4;
  if (remnant > 0)
  {
    encoded[0] = uint8_packing(sequence.substr(0, remnant));
    encoded += 1;
  }

  size_t nb_uint_needed = size / 4;
  for (size_t i = 0; i < nb_uint_needed; i++)
  {
    encoded[i] = uint8_packing(sequence.substr(remnant + 4 * i, 4));
  }
}

using kff_t = std::shared_ptr<Kff_file>;
using kff_raw_t = std::shared_ptr<Section_Raw>;

inline kff_t get_kff(const std::string& path)
{
  uint8_t encoding[] = {0, 1, 3, 2};
  kff_t f = std::make_shared<Kff_file>(path, "w");
  f->write_encoding(encoding);
  return f;
}

inline void write_section(kff_t kff, size_t kmer_size)
{
  Section_GV sgv(kff.get());
  sgv.write_var("k", kmer_size);
  sgv.write_var("max", 1);
  sgv.write_var("data_size", 0);
  sgv.close();
}

inline kff_raw_t get_raw_section(kff_t kff)
{
  return std::make_shared<Section_Raw>(kff.get());
}

template<size_t MAX_K>
inline void write_kmer(kff_raw_t s, KmerSign<MAX_K>& kmer, size_t kmer_size)
{
  std::string k_str = kmer.m_kmer.to_string();
  static uint8_t encoded[1024] = {0};
  encode_sequence(k_str.c_str(), encoded);
  s->write_compacted_sequence(encoded, kmer_size, nullptr);
}

}; // end of namespace kmdiff