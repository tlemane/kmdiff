#pragma once

// std
#include <string>
#include <vector>
#include <fstream>
#include <optional>

// ext
#include <bcli/bcli.hpp>

// int
#include <kmdiff/exceptions.hpp>

namespace kmdiff {

enum class BED_TYPE : uint8_t
{
  BED3 = 3,
  BED4 = 4,
  BED5 = 5,
  BED6 = 6,
  BED12 = 12,
  UNSUPPORTED = 0,
  VISOR_BED = 1
};

struct bed_header
{
  size_t nb_entries;
};

struct bed_record
{
  std::string chrom;
  size_t chrom_start;
  size_t chrom_end;
  std::string name;
  uint16_t score;
  char strand;
  size_t thick_start;
  size_t thick_end;
  std::string item_rgb;
  size_t block_count;
  std::string block_sizes;
  std::string block_starts;
};
typedef struct bed_record bed_record_t;

class BEDReader
{
public:
  BEDReader(const std::string& path);
  void read_all();
  bed_record_t* next();

  auto begin() {return m_records.begin();}
  auto end() {return m_records.end();}
  auto begin() const {return m_records.begin();};
  auto end() const {return m_records.end();}

private:
  void set_bed_type(size_t size);
  void set(const std::string& str_record);

private:
  std::string m_path;
  std::ifstream m_file;
  BED_TYPE m_bed_type;
  bed_record_t m_record;
  std::vector<bed_record_t> m_records;
};

struct bed_record_visor
{
  std::string chrom;
  size_t chrom_start;
  size_t chrom_end;
  std::string type;
  std::string info;
  uint16_t extra;
};
typedef struct bed_record_visor bed_record_visor_t;

class BEDVisorReader
{
public:
  BEDVisorReader(const std::string& path);
  void read_all();
  bed_record_visor_t* next();

  auto begin() {return m_records.begin();}
  auto end() {return m_records.end();}
  auto begin() const {return m_records.begin();};
  auto end() const {return m_records.end();}

private:
  void set_bed_type();

private:
  std::string m_path;
  std::ifstream m_file;
  BED_TYPE m_bed_type;
  bed_record_visor_t m_record;
  std::vector<bed_record_visor_t> m_records;
};

}; // end of namespace kmdiff