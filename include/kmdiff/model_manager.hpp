#pragma once

#include <string>
#include <dlfcn.h>
#include <filesystem>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#define KMTRICKS_PUBLIC
#include <kmtricks/utils.hpp>

#include <kmdiff/exceptions.hpp>

namespace fs = std::filesystem;

namespace kmdiff {

  template<typename Plugin>
  class plugin_manager
  {
    using plugin_t = std::shared_ptr<Plugin>;
    using name_sign_t = std::string (*)();
    using load_sign_t = Plugin* (*)();
    using err_t = const char*;

    private:
      plugin_manager() {}

    public:
      static plugin_manager<Plugin>& get() { static plugin_manager<Plugin> pm; return pm; }

      void init(const std::string& shared_lib_path, const std::string& config)
      {
        m_path = shared_lib_path; m_config = config;

        m_handle = dlopen(m_path.c_str(), RTLD_LAZY);

        if (!m_handle)
          handle_dlerror();

        load_name();
        load_create();
        m_enable = true;

        spdlog::info("Plugin '{}' loaded.", m_name);
      }

      void close()
      {
        if (m_handle)
        {
          dlclose(m_handle);
          m_handle = nullptr;
        }
      }

      plugin_t get_plugin()
      {
        Plugin* p = m_create();
        p->configure(m_config);
        return std::shared_ptr<Plugin>(p);
      }

      std::string name() const
      {
        if (m_enable)
          return m_name;
        return "Not loaded!";
      }

    private:

      void handle_dlerror()
      {
        err_t dlsym_err = dlerror();
        if (dlsym_err)
          throw PluginError(fmt::format("dlerror: {}", dlsym_err));
      }

      void load_create()
      {
        std::string s = std::to_string(sizeof(typename km::selectC<DMAX_C>::type) * 8);
        m_create = reinterpret_cast<load_sign_t>(dlsym(m_handle, fmt::format("create{}", s).c_str()));
        handle_dlerror();
      }

      void load_name()
      {
        name_sign_t plugin_name {nullptr};
        plugin_name = reinterpret_cast<name_sign_t>(dlsym(m_handle, "plugin_name"));
        handle_dlerror();
        m_name = plugin_name();
      }

    private:
      bool m_enable {false};
      bool m_ready {false};
      std::string m_config;
      std::string m_path;
      std::string m_name;

      void* m_handle {nullptr};
      load_sign_t m_create {nullptr};
  };

} // end of namespace kmdiff

