// std
#include <fstream>
#include <string>

// ext
#define _KM_LIB_INCLUDE_
#include <kmtricks/io.hpp>
#include <kmtricks/lz4_stream.hpp>

// int
#include <kmsv/kmer.hpp>

namespace kmsv
{
struct kmer_sign_file_header
{
  uint64_t first_magic;
  uint64_t kmer_size;
  uint64_t partition_id;
  uint64_t is_compressed;
  uint64_t second_magic;
};

typedef struct kmer_sign_file_header kmsheader_t;

const static kmsheader_t kmsheader_d = {km::magic1, 0, 0, km::magic2};

template <typename stream, size_t buf_size = 4096>
class KmerSignFile : public km::IFile<kmsheader_t stream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
  using icstream = lz4_stream::basic_istream<buf_size>;

 public:
  KmerSignFile() = delete;
  KmerSignFile(KmerSignFile const&) = delete;
  KmerSignFile& operator=(KmerSignFile const&) = delete;
  KmerSignFile(KmerSignFile&&) = delete;
  KmerSignFile& operator=(KmerSignFile&&) = delete;

  template <typename S = stream, typename = typename std::enable_if_t<std::is_same<S, km::is>{}, S>>
  KmerSignFile(const std::string& path)
      : km::IFile<kmsheader_t, stream, buf_size>(path, kmsheader_d, std::ios::binary | std::ios::in)
  {
    this->first_layer->read(&this->header, sizeof(this->header));
    if (this->header.first_magic != km::magic1 || this->header.second_magic != km::magic2)
      throw std::runtime_error("Unable to read " + this->path + ". Possibly due to bad format");
    this->template set_second_layer<icstream>(this->header.is_compressed);
  }

  template <typename S = stream, typename = typename std::enable_if_t<std::is_same<S, km::os>{}, S>>
  KmerSignFile(
      const std::string& path, uint64_t kmer_size, uint64_t partition_id, uint64_t compressed)
      : km::IFile<kmsheader_t, stream, buf_size>(
            path, kmsheader_d, std::ios::binary | std::ios::out)
  {
    this->header.kmer_size = kmer_size;
    this->header.partition_id = partition_id;
    this->header.is_compressed = compressed;
    this->first_layer->write(&this->header, sizeof(this->header));
    this->first_layer.flush();
    this->template set_second_layer<ocstream>(this->header.is_compressed);
  }

  template <typename S = stream, typename = typename std::enable_if_t<std::is_same<S, km::os>{}, S>>
  void write(KmerS)
};

};  // namespace kmsv