file(GLOB_RECURSE PLUGIN_SOURCES "*.cpp")

foreach(PLUGIN ${PLUGIN_SOURCES})
  get_filename_component(PLUGIN_LIB ${PLUGIN} NAME_WLE)
  add_library(${PLUGIN_LIB} SHARED ${PLUGIN})
  target_link_libraries(${PLUGIN_LIB} headers)
endforeach()
