cmake_minimum_required(VERSION 3.12)
project(inforf LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(LOCALDEPS_DIR ${PROJECT_BINARY_DIR}/localdeps)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g")


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
