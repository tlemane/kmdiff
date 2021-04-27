#include <gtest/gtest.h>
#include <kmdiff/io.hpp>

using namespace kmdiff;

TEST(vcf, vcf)
{
  {
    VCFReader reader("./data_test/test.vcf");
    bcf1_t* record = reader.next();
    EXPECT_EQ(record->pos, 99);
    EXPECT_EQ(std::string(record->d.allele[0]), "A");
    EXPECT_EQ(std::string(record->d.allele[1]), "T");
    EXPECT_EQ(reader.next(), nullptr);
  }
  {
    VCFReader reader("./data_test/test.vcf");
    reader.read_all();
    int i = 0;
    for (auto record : reader)
    {
      i++;
    }
    EXPECT_EQ(i, 1);
  }

  EXPECT_THROW(VCFReader("./data_test/dummy.txt"), VCFOpenError);
}

TEST(bed, bed)
{
  EXPECT_NO_THROW(BEDReader("./data_test/test3.bed"));
  EXPECT_NO_THROW(BEDReader("./data_test/test4.bed"));
  EXPECT_NO_THROW(BEDReader("./data_test/test5.bed"));
  EXPECT_NO_THROW(BEDReader("./data_test/test6.bed"));
  EXPECT_NO_THROW(BEDReader("./data_test/test12.bed"));
  EXPECT_THROW(BEDReader("./data_test/test9.bed"), BEDBadFormat);
  EXPECT_THROW(BEDReader("./data_test/dummy.txt"), BEDOpenError);
  {
    BEDReader reader("./data_test/test6.bed");
    bed_record_t* record = reader.next();
    EXPECT_EQ(record->name, "A1");
    record = reader.next();
    EXPECT_EQ(record->name, "A2");
    EXPECT_EQ(reader.next(), nullptr);
  }
  {
    BEDReader reader("./data_test/test6.bed");
    reader.read_all();
    int i = 0;
    for (auto record : reader)
    {
      i++;
    }
    EXPECT_EQ(i, 2);
  }
  {
    BEDReader reader("./data_test/test12.bed");
    bed_record_t* record = reader.next();
    EXPECT_EQ(record->chrom, "X");
    EXPECT_EQ(record->chrom_start, 1000);
    EXPECT_EQ(record->chrom_end, 1100);
    EXPECT_EQ(record->name, "A1");
    EXPECT_EQ(record->score, 500);
    EXPECT_EQ(record->strand, '-');
    EXPECT_EQ(record->thick_start, 1020);
    EXPECT_EQ(record->thick_end, 1025);
    EXPECT_EQ(record->item_rgb, "255,0,0");
    EXPECT_EQ(record->block_count, 2);
    EXPECT_EQ(record->block_sizes, "10,10");
    EXPECT_EQ(record->block_starts, "1,1");
  }
}

TEST(bam, bam)
{
  {
    BAMReader reader("./data_test/test.sam");
    bam1_t* record = reader.next();
    EXPECT_EQ(record->core.pos, 41);
    EXPECT_EQ(reader.next(), nullptr);
    EXPECT_EQ(reader.get_header()->n_targets, 1);
  }
  {
    BAMReader reader("./data_test/test.sam");
    reader.read_all();
    int i = 0;
    for (auto record : reader)
    {
      i++;
    }
    EXPECT_EQ(i, 1);
  }
  EXPECT_THROW(BAMReader("./data_test/dummy.txt"), BAMOpenError);
}
