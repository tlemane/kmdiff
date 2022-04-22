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

#pragma once

#include <iterator>
#include <optional>
#include <type_traits>

#include <robin_hood.h>

#include <kmdiff/kmer.hpp>
#include <kmdiff/utils.hpp>

#include <kmtricks/io/lz4_stream.hpp>

#include <spdlog/spdlog.h>

namespace kmdiff {

  template <typename T>
  class IAccumulator
  {
   public:
    IAccumulator() {}

    virtual ~IAccumulator()
    {
    }

    virtual void push(T&& e) = 0;
    virtual size_t size() const = 0;
    virtual void finish() = 0;
    virtual std::optional<T>& get() = 0;
    virtual void destroy() = 0;

   protected:
    std::optional<T> m_opt;
  };

  template <typename T>
  using acc_t = std::shared_ptr<IAccumulator<T>>;

  template <typename T>
  class VectorAccumulator : public IAccumulator<T>
  {
    using iterator = typename std::vector<T>::iterator;

   public:
    VectorAccumulator(size_t reserve = 65536) { m_data.reserve(reserve); }

    void push(T&& e) override
    {
      m_data.push_back(std::move(e));
    }

    void finish() override
    {
      m_it = m_data.begin();
    }

    std::optional<T>& get() override
    {
      if (m_it < m_data.end())
      {
        this->m_opt = *m_it;
        m_it++;
      }
      else
      {
        this->m_opt = std::nullopt;
      }
      return this->m_opt;
    }

    size_t size() const override
    {
      return m_data.size();
    }

    void destroy() override
    {
      std::vector<T>().swap(m_data);
    }

   public:
    std::vector<T> m_data;
    iterator m_it;
  };

  template <typename T>
  class SetAccumulator : public IAccumulator<T>
  {
    using iterator = typename robin_hood::unordered_set<T>::iterator;

   public:
    SetAccumulator(size_t reserve = 65536)
    {
      m_data.reserve(reserve);
    }

    void push(T&& e) override
    {
      m_data.insert(std::move(e));
    }

    void finish() override
    {
      m_it = m_data.begin();
    }

    std::optional<T>& get() override
    {
      if (m_it == m_data.end())
      {
        this->m_opt = std::nullopt;
      }
      else
      {
        this->m_opt = *m_it;
        m_it++;
      }
      return this->m_opt;
    }

    size_t size() const override
    {
      return m_data.size();
    }

    void destroy() override
    {
      robin_hood::unordered_set<T>().swap(m_data);
    }

   private:
    robin_hood::unordered_set<T> m_data;
    iterator m_it;
  };

  template <typename T>
  class FileAccumulator : public IAccumulator<T>
  {
    using o_stream_t = std::ostream;
    using i_stream_t = std::istream;

    using out_stream_t = std::ofstream;
    using in_stream_t = std::ifstream;

    using cpr_out_stream_t = lz4_stream::basic_ostream<8192>;
    using cpr_in_stream_t = lz4_stream::basic_istream<8192>;

   public:
    FileAccumulator(const std::string& path, size_t k_size = 0, bool read = false, bool del = false)
      : m_path(path), m_kmer_size(k_size), m_reading(read), m_del(del)
    {
      if (!m_reading)
      {
        m_out_stream = std::make_shared<out_stream_t>(m_path);
        m_cout_stream = std::make_shared<cpr_out_stream_t>(*m_out_stream);
      }
      else
      {
        m_in_stream = std::make_shared<in_stream_t>(m_path);
        m_cin_stream = std::make_shared<cpr_in_stream_t>(*m_in_stream);
      }

      if constexpr(is_same_template_v<T, KmerSign<32>>)
        m_tmp.set_k(m_kmer_size);
    }

    void push(T&& e) override
    {
      m_size++;
      if constexpr (has_dump<T>::value)
      {
        e.dump(m_cout_stream);
      }
       else
      {
        m_tmp = e;
        m_cout_stream->write(reinterpret_cast<char*>(&m_tmp), sizeof(e));
      }
    }

    void finish() override
    {
      m_cout_stream.reset();
      m_out_stream->close();
      m_out_stream.reset();
      m_in_stream = std::make_shared<in_stream_t>(m_path);
      m_cin_stream = std::make_shared<cpr_in_stream_t>(*m_in_stream);
    }

    std::optional<T>& get() override
    {
      //if (m_read == m_size)
      //{
      //  this->m_opt = std::nullopt;
      //  return this->m_opt;
      //}

      if constexpr (has_load<T>::value)
      {
        bool b = m_tmp.load(m_cin_stream, m_kmer_size);

        if (__builtin_expect(!b, 0))
        {
          this->m_opt = std::nullopt;
          return this->m_opt;
        }

        this->m_opt = m_tmp;
      }
      else
      {
        m_cin_stream->read(reinterpret_cast<char*>(&m_tmp), sizeof(m_tmp));

        if (__builtin_expect(!m_cin_stream->gcount(), 0))
        {
          this->m_opt = std::nullopt;
          return this->m_opt;
        }

        this->m_opt = m_tmp;
      }
      m_read++;
      return this->m_opt;
    }

    size_t size() const override
    {
      return m_size;
    }

    void destroy() override
    {
      if (m_cin_stream)
        m_cin_stream.reset();

      if (m_in_stream)
      {
        m_in_stream->close();
        m_in_stream.reset();
      }

      if (m_del)
      {
        fs::remove(m_path);
      }
    }

    ~FileAccumulator() override
    {
      destroy();
    }

   private:
    T m_tmp;
    size_t m_size{0};
    size_t m_read{0};
    std::string m_path;
    size_t m_kmer_size{0};
    bool m_reading{false};
    bool m_del{false};
    std::shared_ptr<out_stream_t> m_out_stream{nullptr};
    std::shared_ptr<in_stream_t> m_in_stream{nullptr};
    std::shared_ptr<cpr_out_stream_t> m_cout_stream{nullptr};
    std::shared_ptr<cpr_in_stream_t> m_cin_stream{nullptr};
  };

  bool partitions_exist(const std::string& prefix, size_t nb_partitions, const std::string& dir);

}  // end of namespace kmdiff

