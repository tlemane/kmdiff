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

#include <kmdiff/simulator.hpp>

namespace kmdiff
{
std::string convert_to_visor(const std::string& types)
{
  std::vector<std::string> visor_types;
  for (auto& t : bc::utils::split(types, ':'))
  {
    visor_types.push_back(visor_map.at(t));
  }
  return bc::utils::join(visor_types, ",");
}

size_t get_dim_file(const std::string& path, const std::string& output_path)
{
  klibpp::KSeq record;
  klibpp::SeqStreamIn fasta(path.c_str());
  std::ofstream dim(output_path, std::ios::out);
  size_t genome_size = 0;
  if (!dim.good()) throw std::runtime_error(fmt::format("Unable to write at {}", output_path));

  while (fasta >> record)
  {
    dim << record.name << '\t' << std::to_string(record.seq.size()) << "\n";
    genome_size += record.seq.size();
  }
  return genome_size;
}

void get_random_region(
    const std::string& script,
    const std::string& dim_file,
    const std::string& svtype,
    const std::string& svratio,
    const std::string& output_file,
    int n,
    int l,
    int s)
{
  std::string cmd = fmt::format(
      "{} -d {} -n {} -l {} -s {} -v {} -r {} > {}", script, dim_file, n, l, s, svtype, svratio,
      output_file);
  spdlog::debug(cmd);
  std::system(cmd.c_str());
}

void get_mutated(
    const std::string& bin_path,
    const std::string& ref,
    const std::string& bed,
    const std::string& out)
{
  std::string cmd =
      fmt::format("{} HACk -g {} -b {} -o {} > /dev/null 2>&1", bin_path, ref, bed, out);
  spdlog::debug(cmd);
  std::system(cmd.c_str());
}

void get_reads(
    const std::string& ref,
    const std::string& out_template,
    const std::string& snp_out,
    size_t reference_size,
    size_t read_size,
    size_t coverage,
    double error_rate,
    double mutation_rate,
    double indel_fraction,
    double extend,
    size_t seed)
{
  size_t nb_reads = coverage * reference_size / read_size / 2;
  std::string read1 = fmt::format(out_template, 1, seed);
  std::string read2 = fmt::format(out_template, 2, seed);

  Timer t;
  spdlog::debug("Generate reads for {}...", read1);

  FILE* r1 = fopen(read1.c_str(), "w");
  FILE* r2 = fopen(read2.c_str(), "w");
  FILE* snp = fopen(snp_out.c_str(), "w");

  wgsim_core(
      r1, r2, snp, ref.c_str(), 0, nb_reads, 500, 50, read_size, read_size, seed, error_rate,
      mutation_rate, indel_fraction, extend, 0.05);

  fclose(r1);
  fclose(r2);
  fclose(snp);

  spdlog::debug("Generated -> {} ({} seconds)", read1, t.elapsed<std::chrono::seconds>().count());
}

SVPool::SVPool(const std::string& bed_path) : m_bed_path(bed_path) { load(); }

void SVPool::load()
{
  std::ifstream bed(m_bed_path, std::ios::in);
  if (!bed.good()) throw std::runtime_error(fmt::format("Unable to read {}.", m_bed_path));

  for (std::string line; std::getline(bed, line);)
  {
    std::vector<std::string> entries = bc::utils::split(line, '\t');
    m_svs.push_back(
        SV(entries[0], entries[3], entries[4], bc::utils::lexical_cast<size_t>(entries[1]),
           bc::utils::lexical_cast<size_t>(entries[2])));
    m_size++;
  }
}

std::vector<SV> SVPool::get(size_t size, std::mt19937 g)
{
  if (!size) return std::vector<SV>();
  std::shuffle(m_svs.begin(), m_svs.end(), g);
  return slice(m_svs, 0, size + 1);
}

Simulator::Simulator(
    SVPool& control_pool,
    size_t nb_control,
    size_t mean_control,
    size_t sd_control,
    double prob_case,
    SVPool& case_pool,
    size_t nb_case,
    size_t mean_case,
    size_t sd_case,
    double prob_control)
    : m_control_pool(control_pool),
      m_nb_control(nb_control),
      m_mean_control(mean_control),
      m_sd_control(sd_control),
      m_prob_case(prob_case),
      m_case_pool(case_pool),
      m_nb_case(nb_case),
      m_mean_case(mean_case),
      m_sd_case(sd_case),
      m_prob_control(prob_control)
{
}

void Simulator::generate(const std::string& control_dir, const std::string& case_dir)
{
  m_norm_dist = std::normal_distribution<double>(
      static_cast<double>(m_mean_control), static_cast<double>(m_sd_control));
  for (size_t i = 0; i < m_nb_control; i++)
  {
    size_t nb_sv = std::floor(m_norm_dist(m_g));

    spdlog::debug("control{}, nb_sv={}, prob={}", i, nb_sv, m_prob_case);
    auto [nb_control, nb_case] = sample(nb_sv, m_prob_case);
    spdlog::debug("control{}, good={}, bad={}", i, nb_control, nb_case);

    std::string out = fmt::format("{}/control{}", control_dir, i);

    fs::create_directory(out);

    std::vector<SV> from_control = m_control_pool.get(nb_control, m_g);
    std::vector<SV> from_case = m_case_pool.get(nb_case, m_g);

    for (auto& sv : from_case) m_real_control_pool.insert(sv);
    for (auto& sv : from_control) m_real_control_pool.insert(sv);

    spdlog::debug(
        "control{}, from_control_size={}, from_case_size={}", i, from_control.size(),
        from_case.size());

    std::string bed_path = fmt::format("{}/control{}.bed", out, i);
    std::ofstream out_bed(bed_path);
    if (!out_bed.good()) throw std::runtime_error(fmt::format("Unable to write at {}.", bed_path));

    for (auto& sv : from_control) out_bed << sv.to_bed_entry() << "\n";
    for (auto& sv : from_case) out_bed << sv.to_bed_entry() << "\n";
  }

  for (size_t i = 0; i < m_nb_case; i++)
  {
    size_t nb_sv = std::floor(m_norm_dist(m_g));

    spdlog::debug("case{}, nb_sv={}, prob={}", i, nb_sv, m_prob_control);
    auto [nb_case, nb_control] = sample(nb_sv, m_prob_control);
    spdlog::debug("case{}, good={}, bad={}", i, nb_case, nb_control);

    std::string out = fmt::format("{}/case{}", case_dir, i);

    fs::create_directory(out);

    std::vector<SV> from_case = m_case_pool.get(nb_case, m_g);
    std::vector<SV> from_control = m_control_pool.get(nb_control, m_g);

    for (auto& sv : from_case) m_real_case_pool.insert(sv);
    for (auto& sv : from_control) m_real_case_pool.insert(sv);

    spdlog::debug(
        "case{}, from_case_size={}, from_control_size={}", i, from_case.size(),
        from_control.size());

    std::string bed_path = fmt::format("{}/case{}.bed", out, i);
    std::ofstream out_bed(bed_path);
    if (!out_bed.good()) throw std::runtime_error(fmt::format("Unable to write at {}.", bed_path));

    for (auto& sv : from_control) out_bed << sv.to_bed_entry() << "\n";
    for (auto& sv : from_case) out_bed << sv.to_bed_entry() << "\n";
  }

  std::set_intersection(
      m_real_control_pool.begin(), m_real_control_pool.end(), m_real_case_pool.begin(),
      m_real_case_pool.end(), std::back_inserter(m_shared_pool));

  spdlog::debug("control real pool size = {}", m_real_control_pool.size());
  spdlog::debug("case real pool size = {}", m_real_case_pool.size());
  spdlog::debug("nb shared = {}", m_shared_pool.size());
}

void Simulator::generate_refs(
    const std::string& visor_bin,
    const std::string& ref,
    const std::string& control_dir,
    const std::string& case_dir)
{
  for (size_t i = 0; i < m_nb_control; i++)
  {
    std::string path = fmt::format("{}/control{}/control{}.bed", control_dir, i, i);
    std::string out = fmt::format("{}/control{}/control{}_visor", control_dir, i, i);
    get_mutated(visor_bin, ref, path, out);
  }
  for (size_t i = 0; i < m_nb_case; i++)
  {
    std::string path = fmt::format("{}/case{}/case{}.bed", case_dir, i, i);
    std::string out = fmt::format("{}/case{}/case{}_visor", case_dir, i, i);
    get_mutated(visor_bin, ref, path, out);
  }
}

std::string Simulator::generate_reads(
    const std::string& control_dir,
    const std::string& case_dir,
    size_t genome_size,
    size_t read_size,
    size_t coverage,
    double error_rate,
    double mutation_rate,
    double indel_fraction,
    double extend,
    size_t threads)
{
  std::stringstream ss;
  ThreadPool pool(threads);

  for (size_t i = 0; i < m_nb_control; i++)
  {
    std::string ref = fmt::format("{}/control{}/control{}_visor/h1.fa", control_dir, i, i);
    std::string read_dir = fmt::format("{}/control{}/reads", control_dir, i);
    fs::create_directory(read_dir);
    std::string out = read_dir + "/" + "control" + std::to_string(i) + "_{}_{}.fastq";
    std::string snp = fmt::format("{}/control{}_wgsim.bed", read_dir, i);
    size_t seed = m_dist_int(m_g);
    ss << "control" << std::to_string(i) << " : ";
    ss << fmt::format(out, 1, seed) << " ; " << fmt::format(out, 2, seed) << "\n";

    auto task = [=](int id) {
      get_reads(
          ref, out, snp, genome_size, read_size, coverage, error_rate, mutation_rate,
          indel_fraction, extend, seed);
    };
    pool.add_task(task);
  }

  for (size_t i = 0; i < m_nb_case; i++)
  {
    std::string ref = fmt::format("{}/case{}/case{}_visor/h1.fa", case_dir, i, i);
    std::string read_dir = fmt::format("{}/case{}/reads", case_dir, i);
    fs::create_directory(read_dir);
    std::string out = read_dir + "/" + "case" + std::to_string(i) + "_{}_{}.fastq";
    std::string snp = fmt::format("{}/case{}_wgsim.bed", read_dir, i);
    size_t seed = m_dist_int(m_g);
    ss << "case" << std::to_string(i) << " : ";
    ss << fmt::format(out, 1, seed) << " ; " << fmt::format(out, 2, seed) << "\n";
    auto task = [=](int id) {
      get_reads(
          ref, out, snp, genome_size, read_size, coverage, error_rate, mutation_rate,
          indel_fraction, extend, seed);
    };
    pool.add_task(task);
  }
  pool.join_all();
  return ss.str();
}

std::tuple<size_t, size_t> Simulator::sample(size_t size, double prob_bad)
{
  size_t good = 0, bad = 0;
  for (size_t i = 0; i < size; i++)
  {
    double p = m_dist(m_g);
    if (p < prob_bad)
      bad++;
    else
      good++;
  }
  return std::make_tuple(good, bad);
}

void Simulator::dump(
    const std::string& real_control, const std::string& real_case, const std::string& shared)
{
  std::ofstream control_out(real_control, std::ios::out);
  if (!control_out.good())
    throw std::runtime_error(fmt::format("Unable to write at {}.", real_control));
  std::ofstream case_out(real_case, std::ios::out);
  if (!case_out.good()) throw std::runtime_error(fmt::format("Unable to write at {}.", real_case));
  std::ofstream shared_out(shared, std::ios::out);
  if (!shared_out.good()) throw std::runtime_error(fmt::format("Unable to write at {}.", shared));

  for (auto& sv : m_real_control_pool) control_out << sv.to_bed_entry() << "\n";
  for (auto& sv : m_real_case_pool) case_out << sv.to_bed_entry() << "\n";

  for (auto& sv : m_shared_pool) shared_out << sv.to_bed_entry() << "\n";
}

};  // end of namespace kmdiff