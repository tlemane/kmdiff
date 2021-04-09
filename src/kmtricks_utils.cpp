#include <kmdiff/kmtricks_utils.hpp>

namespace kmdiff
{
kmtricks_config_t get_kmtricks_config(const std::string& run_dir)
{
  kmtricks_config_t config;
  const std::string config_path = fmt::format("{}/config.log", run_dir);
  if (!fs::exists(config_path))
    throw KmtricksFileNotFound(fmt::format("kmtricks config file is missing: {}", config_path));
  std::ifstream in_config(config_path, std::ios::in);
  for (std::string line; std::getline(in_config, line);)
  {
    if (bc::utils::contains(line, "Kmer size"))
      config.kmer_size = bc::utils::lexical_cast<size_t>(bc::utils::split(line, ':')[1]);
    if (bc::utils::contains(line, "Nb partitions"))
      config.nb_partitions = bc::utils::lexical_cast<size_t>(bc::utils::split(line, ':')[1]);
  }

  if (!config.kmer_size || !config.nb_partitions)
    throw ConfigError(fmt::format("Unable to load config from {}.", config_path));

  return config;
}

std::vector<std::string> get_fofs(const std::string& run_dir)
{
  std::vector<std::string> fofs;
  kmtricks_config_t config = get_kmtricks_config(run_dir);

  fof_t kmtricks_in = parse_km_fof(fmt::format("{}/storage/fof.txt", run_dir));
  for (int i = 0; i < config.nb_partitions; i++)
  {
    std::string fof_part =
        fmt::format("{}/storage/kmers_partitions/partition_{}/partition{}.fof", run_dir, i, i);

    std::ofstream outfp(fof_part, std::ios::out);
    for (auto& entry : kmtricks_in)
    {
      outfp << fmt::format(
          "{}/storage/kmers_partitions/partition_{}/{}.kmer.lz4\n", run_dir, i, std::get<0>(entry));
    }
    fofs.push_back(fof_part);
  }
  return fofs;
}

std::tuple<std::vector<size_t>, std::vector<size_t>> get_total_kmer(
    const std::string& run_dir, size_t nb_controls, size_t nb_cases)
{
  std::vector<size_t> controls;
  std::vector<size_t> cases;
  std::string path = fmt::format("{}/storage/kmers.total", run_dir);
  std::ifstream in(path, std::ios::in);

  std::string line;
  for (size_t i = 0; i < nb_controls; i++)
  {
    if (!in.good())  // TODO remove
      controls.push_back(100);
    else
    {
      std::getline(in, line);
      controls.push_back(bc::utils::lexical_cast<size_t>(bc::utils::split(line, ' ')[1]));
    }
  }
  for (size_t i = 0; i < nb_cases; i++)
  {
    if (!in.good())
      cases.push_back(100);
    else
    {
      std::getline(in, line);
      cases.push_back(bc::utils::lexical_cast<size_t>(bc::utils::split(line, ' ')[1]));
    }
  }
  return std::make_tuple(controls, cases);
}

};  // namespace kmdiff