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

// std
#include <iterator>
#include <optional>
#include <type_traits>

// ext
#include <robin_hood.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/lz4_stream.hpp>

// int
#include <kmdiff/kmer.hpp>
#include <kmdiff/utils.hpp>

namespace kmdiff
{
template <typename T>
class IAccumulator
{
 public:
  IAccumulator() {}

  virtual void push(T&& e) = 0;
  virtual size_t size() const = 0;
  virtual void finish() = 0;
  virtual std::optional<T>& get() = 0;

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

  void push(T&& e) override { m_data.push_back(std::move(e)); }

  void finish() override { m_it = m_data.begin(); }

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

  size_t size() const override { return m_data.size(); }

 public:
  std::vector<T> m_data;
  iterator m_it;
};

template <typename T>
class SetAccumulator : public IAccumulator<T>
{
  using iterator = typename robin_hood::unordered_set<T>::iterator;

 public:
  SetAccumulator(size_t reserve = 65536) { m_data.reserve(reserve); }

  void push(T&& e) override { m_data.insert(std::move(e)); }

  void finish() override { m_it = m_data.begin(); }

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

  size_t size() const override { return m_data.size(); }

 private:
  robin_hood::unordered_set<T> m_data;
  iterator m_it;
};

template <typename T>
class FileAccumulator : public IAccumulator<T>
{
  // using out_stream_t = lz4_stream::basic_ostream<8192>;
  // using in_stream_t = lz4_stream::basic_istream<8192>;
  using out_stream_t = std::ofstream;
  using in_stream_t = std::ifstream;

 public:
  FileAccumulator(const std::string& path, size_t k_size = 0) : m_path(path), m_kmer_size(k_size)
  {
    m_out_stream = std::make_shared<out_stream_t>(m_path);
  }

  void push(T&& e) override
  {
    m_size++;
    if constexpr (has_dump<T>::value)
      e.dump(m_out_stream);
    else
    {
      m_tmp = e;
      m_out_stream->write(reinterpret_cast<char*>(&m_tmp), sizeof(e));
    }
  }

  void finish() override
  {
    m_out_stream->close();
    m_out_stream.reset();
    m_in_stream = std::make_shared<in_stream_t>(m_path);
  }

  std::optional<T>& get() override
  {
    if (m_read == m_size)
    {
      this->m_opt = std::nullopt;
      return this->m_opt;
    }

    if constexpr (has_load<T>::value)
    {
      m_tmp.load(m_in_stream, m_kmer_size);
      this->m_opt = m_tmp;
    }
    else
    {
      m_in_stream->read(reinterpret_cast<char*>(&m_tmp), sizeof(m_tmp));
      this->m_opt = m_tmp;
    }
    m_read++;
    return this->m_opt;
  }

  size_t size() const override { return m_size; }

 private:
  T m_tmp;
  size_t m_size{0};
  size_t m_read{0};
  size_t m_kmer_size{0};
  std::string m_path;
  std::shared_ptr<out_stream_t> m_out_stream{nullptr};
  std::shared_ptr<in_stream_t> m_in_stream{nullptr};
};

};  // end of namespace kmdiff