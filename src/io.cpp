#include <kmdiff/io.hpp>

namespace kmdiff {

VCFReader::VCFReader(const std::string& path)
  : m_path(path)
{
  m_record = bcf_init();
  m_vcf_file = bcf_open(m_path.c_str(), "r");
  if (m_vcf_file == NULL)
    throw VCFOpenError(fmt::format("Unable to open {}.", m_path));
  m_header = bcf_hdr_read(m_vcf_file);
  if (m_header == NULL)
    throw VCFHeaderError(fmt::format("Unable to read header: {}", m_path));
}

VCFReader::~VCFReader()
{
  bcf_hdr_destroy(m_header);
  bcf_destroy(m_record);
  bcf_close(m_vcf_file);
  for (auto& r : m_records)
    bcf_destroy(r);
}

void VCFReader::read_all()
{
  while (bcf_read(m_vcf_file, m_header, m_record) == 0)
    m_records.push_back(bcf_dup(m_record));
}

bcf1_t* VCFReader::next()
{
  if (bcf_read(m_vcf_file, m_header, m_record) == 0)
  {
    bcf_unpack(m_record, BCF_UN_ALL);
    return m_record;
  }
  return nullptr;
}

BEDReader::BEDReader(const std::string& path)
  : m_path(path)
{
  m_file = std::ifstream(path, std::ios::in);
  if (!m_file.good())
    throw BEDOpenError(fmt::format("Unable to open {}.", m_path));
  for (std::string line; std::getline(m_file, line);)
  {
    if (!bc::utils::startswith(line, "#"))
    {
      set_bed_type(bc::utils::split(line, '\t').size());
      break;
    }
  }
  if (m_bed_type == BED_TYPE::UNSUPPORTED)
    throw BEDBadFormat(fmt::format("Bed bad format: {}", m_path));
  m_file.seekg(std::ios::beg);
}

void BEDReader::set_bed_type(size_t size)
{
  switch (size)
  {
  case 3:
    m_bed_type = BED_TYPE::BED3; break;
  case 4:
    m_bed_type = BED_TYPE::BED4; break;
  case 5:
    m_bed_type = BED_TYPE::BED5; break;
  case 6:
    m_bed_type = BED_TYPE::BED6; break;
  case 12:
    m_bed_type = BED_TYPE::BED12; break;
  default:
    m_bed_type = BED_TYPE::UNSUPPORTED; break;
  }
}

void BEDReader::set(const std::string& str_record)
{
  std::vector<std::string> s = bc::utils::split(str_record, '\t');
  if (static_cast<uint8_t>(m_bed_type) >= 3)
  {
    m_record.chrom = s[0];
    m_record.chrom_start = bc::utils::lexical_cast<size_t>(s[1]);
    m_record.chrom_end = bc::utils::lexical_cast<size_t>(s[2]);
  }

  if (static_cast<uint8_t>(m_bed_type) >= 4)
    m_record.name = s[3];

  if (static_cast<uint8_t>(m_bed_type) >= 5)
    m_record.score = bc::utils::lexical_cast<uint16_t>(s[4]);

  if (static_cast<uint8_t>(m_bed_type) == 12)
  {
    m_record.strand = bc::utils::lexical_cast<char>(s[5]);
    m_record.thick_start = bc::utils::lexical_cast<size_t>(s[6]);
    m_record.thick_end = bc::utils::lexical_cast<size_t>(s[7]);
    m_record.item_rgb = s[8];
    m_record.block_count = bc::utils::lexical_cast<size_t>(s[9]);
    m_record.block_sizes = s[10];
    m_record.block_starts = s[11];
  }
}

void BEDReader::read_all()
{
  for (std::string line; std::getline(m_file, line);)
  {
    if (!bc::utils::startswith(line, "#"))
    {
      set(line);
      m_records.push_back(m_record);
    }
  }
}

bed_record_t* BEDReader::next()
{
  std::string line;
  if (std::getline(m_file, line))
  {
    set(line);
    return &m_record;
  }
  return nullptr;
}

BAMReader::BAMReader(const std::string& path)
  : m_path(path)
{
  m_bam_file = hts_open(m_path.c_str(), "r");
  if (m_bam_file == NULL)
    throw BAMOpenError(fmt::format("Unable to open {}", m_path));
  m_bam_header = sam_hdr_read(m_bam_file);
  if (m_bam_header == NULL)
    throw BAMHeaderError(fmt::format("Unable to read bam header: {}", m_path));
  m_record = bam_init1();
}

BAMReader::~BAMReader()
{
  bam_hdr_destroy(m_bam_header);
  bam_destroy1(m_record);
  sam_close(m_bam_file);
  for (auto& r : m_records)
    bam_destroy1(r);
}

bam_hdr_t* BAMReader::get_header()
{
  return m_bam_header;
}

void BAMReader::read_all()
{
  while(sam_read1(m_bam_file, m_bam_header, m_record) == 0)
  {
    m_records.push_back(bam_dup1(m_record));
  }
}

bam1_t* BAMReader::next()
{
  int res = sam_read1(m_bam_file, m_bam_header, m_record);
  if (res == 0)
    return m_record;
  return nullptr;
}

}; // end of namespace kmdiff