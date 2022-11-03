/*****************************************************************************
 *   kmdiff
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <cstdint>
#include <string>
#include <memory>

#include <kff_io.hpp>

#include <kmdiff/kmer.hpp>

namespace kmdiff {

  using kff_t = std::unique_ptr<Kff_file>;
  using kff_raw_t = std::unique_ptr<Section_Raw>;

  class KffWriter
  {
  public:
    KffWriter(const std::string& path, size_t kmer_size)
      : m_kmer_size(kmer_size)
    {
      m_kff_file = std::make_unique<Kff_file>(path, "w");
      uint8_t encoding[] = {0, 1, 3, 2};
      m_kff_file->write_encoding(encoding);

      Section_GV sgv(m_kff_file.get());
      sgv.write_var("k", m_kmer_size);
      sgv.write_var("max", 1);
      sgv.write_var("data_size", 0);
      sgv.close();
      m_kff_sec = std::make_unique<Section_Raw>(m_kff_file.get());
    }

    template<size_t MAX_K>
    void write(KmerSign<MAX_K>& kmer)
    {
      uint8_t encoded[1024];
      m_kmer = kmer.m_kmer.to_string();
      encode_sequence(encoded);
      m_kff_sec->write_compacted_sequence(encoded, m_kmer_size, nullptr);
    }

    template<size_t MAX_K>
    void write(km::Kmer<MAX_K>& kmer)
    {
      uint8_t encoded[1024];
      m_kmer = kmer.to_string();
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

  using kff_reader_t = std::unique_ptr<Kff_reader>;

  class KffReader
  {
    std::string nt[4] = {"A", "C", "T", "G"};
  public:
    KffReader(const std::string& path, size_t kmer_size)
    {
      m_kff_reader = std::make_unique<Kff_reader>(path);
      m_kmer_size = kmer_size;
      m_data_size = 0;
      for (int i=0; i<256; i++)
      {
        for (int j=0; j<4; j++)
        {
          uint8_t n = (i>>(2*j)) & 0b11;
          m_lookup[i] = nt[n] + m_lookup[i];
        }
      }
    }

    template<size_t MAX_K>
    std::optional<km::Kmer<MAX_K>> read()
    {
      static km::Kmer<MAX_K> kmer;
      kmer.set_k(m_kmer_size);
      if (m_kff_reader->has_next())
      {
        m_kff_reader->next_kmer(m_buffer, m_data);
        return km::Kmer<MAX_K>(to_string());
      }
      return std::nullopt;
    }

  private:
    std::string to_string()
    {
      size_t size = m_kmer_size % 4 == 0 ? m_kmer_size / 4 : m_kmer_size / 4 + 1;
      std::string str_kmer = m_lookup[m_buffer[0]];
      if (m_kmer_size % 4 != 0)
      {
        int trunc = 4 - (m_kmer_size % 4);
        str_kmer = str_kmer.substr(trunc);
      }
      for (size_t i=1; i<size; i++)
      {
        uint8_t v = m_buffer[i];
        str_kmer += m_lookup[v];
      }
      return str_kmer;
    }

    kff_reader_t m_kff_reader {nullptr};
    size_t m_kmer_size{0};
    size_t m_data_size{0};
    uint8_t* m_buffer {nullptr};
    uint8_t* m_data {nullptr};
    std::string m_lookup[256];
  };

  using kff_r_t = std::unique_ptr<KffReader>;

} // end of namespace kmdiff

