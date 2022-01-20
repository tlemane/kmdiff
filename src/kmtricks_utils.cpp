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

#include <kmdiff/kmtricks_utils.hpp>

#define KMTRICKS_PUBLIC
#include <kmtricks/io/fof.hpp>
#include <kmtricks/kmdir.hpp>
#include <kmtricks/io/hist_file.hpp>

namespace kmdiff {

kmtricks_config_t get_kmtricks_config(const std::string& run_dir)
{
  kmtricks_config_t config;

  const std::string config_path = fmt::format("{}/kmdiff-count.opt", run_dir);

  std::ifstream in_config(config_path, std::ios::in);
  for (std::string line; std::getline(in_config, line);)
  {
    if (bc::utils::contains(line, "kmer_size"))
    {
      for (auto&& o : bc::utils::split(line, ','))
      {
        if (bc::utils::contains(o, "kmer_size"))
          config.kmer_size = bc::utils::lexical_cast<size_t>(bc::utils::split(o, '=')[1]);
        if (bc::utils::contains(o, "abundance_min"))
          config.abundance_min = bc::utils::lexical_cast<size_t>(bc::utils::split(o, '=')[1]);
      }
    }
  }

  std::size_t np = 0;
  for (auto& _ : fs::directory_iterator(fmt::format("{}/counts", run_dir)))
    np++;
  config.nb_partitions = np;

  if (!config.kmer_size || !config.nb_partitions)
    throw ConfigError(fmt::format("Unable to load config from {}.", config_path));

  return config;
}


km::Fof get_fofs(const std::string& run_dir)
{
  std::string fof_path = fmt::format("{}/kmtricks.fof", run_dir);
  return km::Fof(fof_path);
}

std::tuple<std::vector<size_t>, std::vector<size_t>> get_total_kmer(
  const std::string& run_dir,
  size_t nb_controls,
  size_t nb_cases,
  size_t ab_min)
{
  auto fof = get_fofs(run_dir);
  std::vector<size_t> total_controls(nb_controls);
  std::vector<size_t> total_cases(nb_cases);

  for (std::size_t i = 0; i < total_controls.size(); i++)
  {
    std::string fid = fof.get_id(i);
    std::string hpath = km::KmDir::get().get_hist_path(fid);
    km::HistReader hr(hpath);

    auto hist = hr.get();
    auto info = hr.infos();
    auto& v = hist->get_vec();

    total_controls[i] = info.total;
    for (std::size_t j=1; j<ab_min; j++)
      total_controls[i] -= (j) * v[j-1];

    spdlog::debug("{}: {} k-mers", fid, total_controls[i]);
  }

  for (std::size_t i = nb_controls; i < nb_controls + nb_cases; i++)
  {
    std::string fid = fof.get_id(i);
    std::string hpath = km::KmDir::get().get_hist_path(fid);
    km::HistReader hr(hpath);

    auto hist = hr.get();
    auto info = hr.infos();
    auto& v = hist->get_vec();

    total_cases[i - nb_controls] = info.total;
    for (std::size_t j=1; j<ab_min; j++)
      total_cases[i - nb_controls] -= (j) * v[j-1];

    spdlog::debug("{}: {} k-mers", fid, total_cases[i - nb_controls]);
  }

  return std::make_tuple(std::move(total_controls), std::move(total_cases));
}

} // end of namespace kmdiff
