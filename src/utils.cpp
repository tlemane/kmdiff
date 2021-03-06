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

#include <random>

#include <spdlog/spdlog.h>
#include <bcli/bcli.hpp>

#include <kmdiff/utils.hpp>

#include <sys/wait.h>
#include <sys/resource.h>

#if __APPLE__
  #include <mach-o/dyld.h>
  #include <mach/mach.h>
#endif

namespace kmdiff {

  std::string command_exists(const std::string& dir, const std::string& cmd)
  {
    if (!std::system(fmt::format("which {} > /dev/null 2>&1", cmd).c_str()))
      return cmd;
    else if (!std::system(fmt::format("which {}/{} > /dev/null 2>&1", dir, cmd).c_str()))
      return fmt::format("{}/{}", dir, cmd);
    else
      throw BinaryNotFound(fmt::format("{} not found.", cmd));
  }

  std::string get_binary_dir()
  {
    uint32_t size = 1024;
    char buffer[1024];
  #if __APPLE__
    _NSGetExecutablePath(buffer, &size);
  #else
    auto r = readlink("/proc/self/exe", buffer, size);
    if (r == -1)
      throw BinaryNotFound("Unable to found kmdiff binary path.");
  #endif
    return fs::path(buffer).parent_path();
  }

  std::string get_uname_sr()
  {
    std::array<char, 256> buffer;
    std::string result;
    auto pipe = popen("uname -sr", "r");

    while (!feof(pipe))
    {
      if (fgets(buffer.data(), 256, pipe) != nullptr) result += buffer.data();
    }
    auto rc = pclose(pipe); unused(rc);

    return result;
  }

  VerbosityLevel str_to_verbosity_level(const std::string& str_level)
  {
    if (str_level == "debug")
      return VerbosityLevel::DEBUG;
    else if (str_level == "info")
      return VerbosityLevel::INFO;
    else if (str_level == "warning")
      return VerbosityLevel::WARNING;
    else if (str_level == "error")
      return VerbosityLevel::ERROR;
    else
      return VerbosityLevel::WARNING;
  }

  void set_verbosity_level(const std::string& level)
  {
    switch (str_to_verbosity_level(level))
    {
      case VerbosityLevel::DEBUG:
        spdlog::set_level(spdlog::level::debug);
        break;
      case VerbosityLevel::INFO:
        spdlog::set_level(spdlog::level::info);
        break;
      case VerbosityLevel::WARNING:
        spdlog::set_level(spdlog::level::warn);
        break;
      case VerbosityLevel::ERROR:
        spdlog::set_level(spdlog::level::err);
        break;
    }
  }

  int exec_external_cmd(const std::string& cmd, const std::string& args,
                        const std::string& sout, const std::string& serr)
  {
    std::string bin_name = fs::path(cmd).filename().string();
    std::vector<std::string> argsv = bc::utils::split(args, ' ');
    pid_t pid;
    int p[2];
    pipe(p);
    int status;
    const char* cmd_name[] = {cmd.c_str()};
    const char** argv = new const char*[argsv.size() + 2];
    for (size_t i = 1; i < argsv.size() + 1; i++)
    {
      argv[i] = argsv[i - 1].c_str();
    }

    spdlog::info(fmt::format("exec: {} {}", bin_name, args));
    argv[0] = cmd_name[0];
    argv[argsv.size() + 1] = NULL;

    pid = fork();
    if (pid == 0)
    {
      int fd1 = -1, fd2 = -1;

      if (!sout.empty())
      {
        fd1 = open(sout.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
        dup2(fd1, STDOUT_FILENO);
      }

      if (!serr.empty())
      {
        fd2 = open(serr.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
        dup2(fd2, STDERR_FILENO);
      }

      if (bool failed = execvp(cmd_name[0], (char**)argv) == -1; failed)
      {
        delete[] argv;
        write(p[1], &failed, sizeof(failed));
        throw ExternalExecFailed(fmt::format("Failed to run {}.", cmd_name[0]));
      }

      close(fd1);
      close(fd2);

      delete[] argv;
    }
    else if (pid > 0)
    {
      if (wait(&status) == -1)
      {
        throw ExternalExecFailed("wait() failed.");
      }
      fcntl(p[0], F_SETFL, fcntl(p[0], F_GETFL) | O_NONBLOCK);
      bool failed = false;
      ssize_t r = read(p[0], &failed, sizeof(failed)); unused(r);

      if (WIFEXITED(status))
        if (WEXITSTATUS(status) != 0)
          throw ExternalExecFailed(fmt::format("{} exit with {}.", cmd_name[0], WEXITSTATUS(status)));

      if (!failed)
      {
        spdlog::info("{} exit normally with ({}).", bin_name, WEXITSTATUS(status));
        return 0;
      }
      else
        return 1;
    }
    else
    {
      throw ExternalExecFailed("fork() failed.");
    }
    return 1;
  }

  std::string random_dna_seq(size_t size)
  {
    static const char alpha[] = "ACGT";
    static std::default_random_engine g;
    static std::uniform_int_distribution<int> dist(0, 3);
    std::string seq;
    for (size_t i = 0; i < size; i++) seq.push_back(alpha[dist(g)]);
    return seq;
  }

  size_t get_peak_rss()
  {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
  #if __APPLE__
    return static_cast<size_t>(rusage.ru_maxrss);
  #else
    return static_cast<size_t>(rusage.ru_maxrss);
  #endif
  }

  size_t get_current_rss()
  {
  #if __APPLE__
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
      (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
      return (size_t)0L;
    return (size_t)info.resident_size;
  #else
    long rss = 0;
    FILE *fp = NULL;
    if ((fp = fopen("/proc/self/statm", "r")) == NULL)
      return 0;
    if (fscanf(fp, "%*s%ld", &rss) != 1)
    {
      fclose(fp);
      return 0;
    }
    fclose(fp);
    return rss * sysconf(_SC_PAGE_SIZE);
  #endif
  }

  std::string& str_to_upper(std::string& s)
  {
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
  }

} // end of namespace kmdiff

