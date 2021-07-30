cmake_minimum_required(VERSION 3.10)
project(inforf LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(LOCALDEPS_DIR ${PROJECT_BINARY_DIR}/localdeps)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

set(ORTOOLS_ROOT ${ORTOOLS_ROOT} ${CMAKE_SOURCE_DIR}/libs/ortools)

message(${ORTOOLS_ROOT})

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${ORTOOLS_ROOT})

find_package(ortools REQUIRED)

# multithreading for std::async
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


add_subdirectory(libs)
add_subdirectory(src)
add_subdirectory(test)


add_subdirectory(googletest EXCLUDE_FROM_ALL)

include_directories(SYSTEM ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})


#target_link_libraries(rfdist inforf)
#include_directories(SYSTEM ${gtest_SOURCE_DIR}/../googlemock/include ${gtest_SOURCE_DIR}/../googlemock/)

