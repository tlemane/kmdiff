#pragma once

// std
#include <string>
#include <vector>

// ext
#include <htslib/sam.h>
#include <fmt/format.h>

// int
#include <kmdiff/exceptions.hpp>

namespace kmdiff {

class BAMReader
{
public:
  BAMReader(const std::string& path);
  ~BAMReader();

  void read_all();
  bam1_t* next();
  bam_hdr_t* get_header();

  auto begin() {return m_records.begin();}
  auto end() {return m_records.end();}
  auto begin() const {return m_records.begin();};
  auto end() const {return m_records.end();}

private:
  std::string m_path;
  samFile* m_bam_file;
  bam_hdr_t* m_bam_header;
  bam1_t* m_record;
  std::vector<bam1_t*> m_records;
};

}; // end of namespace kmdiff