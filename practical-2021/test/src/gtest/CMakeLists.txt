
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