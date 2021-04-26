file(GLOB test_src
     "*.cpp"
)
add_executable(Google_Tests_run ${test_src})
target_link_libraries(Google_Tests_run
  gtest
  gtest_main
  inforf
  ${PLLMODULES_LIBRARIES}
  ${PLL_LIBRARIES}
  ${PLLMODULES_LIBRARIES}
)

target_include_directories(Google_Tests_run PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${LOCALDEPS_DIR}/include>)

add_test(Google_Tests_run Google_Tests_run)

set(target Google_Tests_run)

add_custom_command(TARGET ${target}
                   POST_BUILD
                   COMMAND ${target}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   COMMENT "Running ${target}" VERBATIM)
