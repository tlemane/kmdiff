file(GLOB CPP_FILES RELATIVE ${PROJECT_SOURCE_DIR}/src "*.cpp")
file(GLOB CPP_FILES_2 RELATIVE ${PROJECT_SOURCE_DIR} "*.cpp")

include(ClangFormat)

set(PROJECT_FILES "")
foreach(project_file ${CPP_FILES_2})
  list(APPEND PROJECT_FILES ${project_file})
endforeach()
foreach(project_file ${HPP_FILES})
  list(APPEND PROJECT_FILES ${project_file})
endforeach()

add_format_project_target("${PROJECT_FILES}")

add_library(${PROJECT_NAME}lib ${CPP_FILES})
target_link_libraries(${PROJECT_NAME}lib PUBLIC build_type_flags headers links deps)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}lib)

if(WITH_POPSTRAT)
  add_compile_definitions(WITH_POPSTRAT)
endif()

if (NOT APPLE)
  target_link_options(${PROJECT_NAME} PUBLIC -export-dynamic)
endif()
if (STATIC_BUILD)
  target_link_options(${PROJECT_NAME} PUBLIC -static)
endif()
