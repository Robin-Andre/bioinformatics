
option(ENABLE_TESTS "Manually Compile Googletest (currently incomapatible with softwipe)" OFF)
if(ENABLE_TESTS) 

  set(TEST_EXE_NAME Google_Test_run)

  message("Testcase compilation active")
  #Collect all testcases into one target (does not support directories yet)
  file(GLOB test_src
     "*.cpp"
  )
  #Compile Google Executable
  add_executable(Google_Test_run ${test_src})
  #This is an attempt to get the git hash into c++ 
  execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest commit hash
  execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  message(Git Hash ${GIT_COMMIT_HASH})
  message("I do not know what to include")
  message(${CMAKE_SOURCE_DIR})
  target_include_directories(${TEST_EXE_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/src/)
  target_compile_definitions(${TEST_EXE_NAME} PRIVATE "-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")
  #googletest requires gtest as library. gtest is not linked as system library which causes a lot of
  #-Wundef errors in softwipe which are 
  target_link_libraries(${TEST_EXE_NAME} gtest_main)
  #gtest now requires gmp 
  target_link_libraries(${TEST_EXE_NAME} gmp :libgmp.a)
  #gtest also requires json_lib
  #target_link_libraries(${TEST_EXE_NAME} json_lib)
  #Linking the entire project might be a bit overkill but the test cases need a link to the classes somewhere
  #maybe figure out a smoother solution to this ... later
  target_link_libraries(${TEST_EXE_NAME} inforf)
  #This includes libpll somehow. I again it might be solved more elegantly but at a later stage
  target_include_directories(${TEST_EXE_NAME} SYSTEM PUBLIC $<BUILD_INTERFACE:${LOCALDEPS_DIR}/include>)  
  #This tells Cmake that this executable is meant as a test (I have not yet understood what it does)
  add_test(Random_Name ${TEST_EXE_NAME})

  #add_executable(t )
  #And this ensures that the tests do run from make
  set(target ${TEST_EXE_NAME})
  add_custom_command(TARGET ${target}
                   POST_BUILD
                   COMMAND ${target}
                   WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
                   COMMENT "Running ${target}" VERBATIM)
else()
  message(WARNING "Testcase compilation not active")
endif()