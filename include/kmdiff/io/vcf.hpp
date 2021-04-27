#pragma once

// std
#include <string>

// ext
#include <htslib/vcf.h>
#include <fmt/format.h>

// int
#include <kmdiff/exceptions.hpp>

namespace kmdiff {

class VCFReader
{
public:
  VCFReader(const std::string& path);
  ~VCFReader();

  void read_all();
  bcf1_t* next();

  auto begin() {return m_records.begin();}
  auto end() {return m_records.end();}
  auto begin() const {return m_records.begin();};
  auto end() const {return m_records.end();};

private:
  std::string m_path;
  htsFile* m_vcf_file;
  bcf_hdr_t* m_header;
  bcf1_t* m_record;
  std::vector<bcf1_t*> m_records;
};

}; // end of namespace kmdiff