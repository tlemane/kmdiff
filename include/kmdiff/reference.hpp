#pragma once

// std
#include <string>

// ext
#include <htslib/faidx.h>
#include <fmt/format.h>

// int
#include <kmdiff/sv.hpp>
#include <kmdiff/kmer.hpp>

namespace kmdiff {

class Reference
{
public:
  Reference (const std::string& path);
  ~Reference();
  std::string fetch(const std::string& name, size_t start, size_t end);
  std::tuple<std::string, std::string> get_kmer_view(SV& sv, size_t kmer_size);

private:
  std::string m_path;
  faidx_t*    m_faidx {nullptr};
};

};