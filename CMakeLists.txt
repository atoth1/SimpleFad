cmake_minimum_required(VERSION 3.18.0 FATAL_ERROR)

project(SimpleFad VERSION 0.0.0)

include(GNUInstallDirs)
include(CTest)
include(GoogleTest)
find_package(GTest REQUIRED)

option(ENABLE_EXCEPTIONS
	"Enable runtime exceptions (bounds checks, valid function inputs, etc.)."
	OFF)

option(ENABLE_TESTING "Enable building of unit tests." ON)

# Tests demonstrating Newton's method utilize LAPACK routines with name mangling hard-coded
# to the function name followed by a single underscore. If LAPACK library symbol names don't
# match this convention these particular tests may be disabled with this option.
option(ENABLE_NEWTON_TESTS
	"Enable building of Newton's method tests which include hard-coded LAPACK symbol names."
	ON)

# Surprised that the gcc standard library on my machine supports this, so why not make
# it an option. Expecially since there look to be some proposals to make more/most of
# <cmath> constexpr in the standard. Doesn't work with clang on my machine.
option(CMATH_SUPPORTS_CONSTEXPR
	"Does <cmath> support (non-standard) constexpr for required mathematical operations."
	OFF)

add_subdirectory(src ${CMAKE_CURRENT_BINARY_DIR}/src)
if (ENABLE_TESTING)
  find_package(LAPACK)
  add_subdirectory(test ${CMAKE_CURRENT_BINARY_DIR}/test)
endif()

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
  
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME})

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}_Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  
install(
  EXPORT ${PROJECT_NAME}_Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME})
  
install(
  FILES
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME})
