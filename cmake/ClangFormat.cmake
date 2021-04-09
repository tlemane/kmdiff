find_program(CLANG_FORMAT_BIN NAMES "clang-format" DOC "Path to clang-format executable") 
if(NOT CLANG_FORMAT_BIN) 
  message(STATUS "clang-format not found.") 
else() 
  message(STATUS "clang-format found: ${CLANG_FORMAT_BIN}. Add target format-project.") 
endif()

macro(add_format_project_target FILES_TO_FORMAT)
  if(CLANG_FORMAT_BIN)
    add_custom_target(
      format-project
      COMMAND ${CLANG_FORMAT_BIN} -i -style=file ${FILES_TO_FORMAT}
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  endif()
endmacro()
