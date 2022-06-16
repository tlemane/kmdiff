SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  ${PROJECT_NAME})
SET (CPACK_PACKAGE_VERSION_MAJOR        "${CMAKE_PROJECT_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${CMAKE_PROJECT_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${CMAKE_PROJECT_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${CMAKE_PROJECT_VERSION}")
SET (CPACK_PACKAGE_FILE_NAME            "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CMAKE_SYSTEM_NAME}")
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")
SET (CPACK_SET_DESTDIR true)
SET (CPACK_INSTALL_PREFIX /)

include (CPack)
