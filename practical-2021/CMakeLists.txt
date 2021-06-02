cmake_minimum_required(VERSION 3.12)
project(inforf LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(LOCALDEPS_DIR ${PROJECT_BINARY_DIR}/localdeps)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

set(ORTOOLS_ROOT ${CMAKE_SOURCE_DIR}/libs/ortools)

message(WARNING ${ORTOOLS_ROOT})

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${ORTOOLS_ROOT})

find_package(ortools REQUIRED)

#Begin Hacky internet locate GMP
include(FeatureSummary)
set_package_properties(GMP PROPERTIES
  DESCRIPTION "GNU multi-precision library"
  URL "https://gmplib.org"
)

# Try finding the package with pkg-config
find_package(PkgConfig QUIET)
pkg_check_modules(PKG QUIET gmp gmpxx)

# Try to locate the libraries and their headers, using pkg-config hints
find_path(GMP_INCLUDE_DIR gmp.h HINTS ${PKG_gmp_INCLUDEDIR})
find_library(GMP_LIB gmp HINTS ${PKG_gmp_LIBDIR})

find_path(GMPXX_INCLUDE_DIR gmpxx.h HINTS ${PKG_gmpxx_INCLUDEDIR})
find_library(GMPXX_LIB gmpxx HINTS ${PKG_gmpxx_LIBDIR})

# Remove these variables from cache inspector
mark_as_advanced(GMP_INCLUDE_DIR GMP_LIB GMPXX_INCLUDE_DIR GMPXX_LIB)

# Report if package was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  DEFAULT_MSG
  GMPXX_LIB GMPXX_INCLUDE_DIR GMP_INCLUDE_DIR GMP_LIB
)

# Set targets
if(GMP_FOUND)
  # C library
  if(NOT TARGET GMP::gmp)
    add_library(GMP::gmp UNKNOWN IMPORTED)
    set_target_properties(GMP::gmp PROPERTIES
      IMPORTED_LOCATION ${GMP_LIB}
      INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE_DIR}
    )
  endif()

  # C++ library, which requires a link to the C library
  if(NOT TARGET GMP::gmpxx)
    add_library(GMP::gmpxx UNKNOWN IMPORTED)
    set_target_properties(GMP::gmpxx PROPERTIES
      IMPORTED_LOCATION ${GMPXX_LIB}
      INTERFACE_INCLUDE_DIRECTORIES ${GMPXX_INCLUDE_DIR}
      INTERFACE_LINK_LIBRARIES GMP::gmp
    )
  endif()
endif()

#End hacky internet locate GMP
#find_package(GMP REQUIRED)

add_subdirectory(libs)
add_subdirectory(src)
add_subdirectory(utils)
add_subdirectory(test)


add_subdirectory(googletest EXCLUDE_FROM_ALL)

include_directories(SYSTEM ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

#target_link_libraries(rfdist inforf)
#include_directories(SYSTEM ${gtest_SOURCE_DIR}/../googlemock/include ${gtest_SOURCE_DIR}/../googlemock/)

#add_custom_target(softwipe ALL
#    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/softwipe.txt")

#set(SOFTWIPE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../softwipe/softwipe_detailed_result/)

#add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/softwipe.txt"
#                 WORKING_DIRECTORY ${SOFTWIPE_DIRECTORY}
#                  COMMAND python3.9 ../softwipe.py -Cm -e ../../practical-2021/softwipe_call -f ../../practical-2021/softwipe_command_file  ../../practical-2021/ -x ../../practical-2021/googletest/,../../practical-2021/libs
#                  COMMENT "Running softwipe"
#                  POST_BUILD)
#python3.9 ../softwipe.py -Cm -e ../../practical-2021/softwipe_call -f ../../practical-2021/softwipe_command_file  ../../practical-2021/ -x ../../practical-2021/googletest/,../../practical-2021/libs,../../practical-2021/utils
#python3.9 ../softwipe.py -Cm -e ../../practical-2021/softwipe_call -f ../../practical-2021/softwipe_command_file ../../practical-2021/src/
