#include <cstdint>
#include <string>

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