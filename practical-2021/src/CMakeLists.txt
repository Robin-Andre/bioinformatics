add_library(inforf STATIC
    datastructures/PllTree.cpp
    datastructures/PllTree.hpp
    datastructures/PllSplits.cpp
    datastructures/PllSplits.hpp
)


set_target_properties(inforf PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

target_link_libraries(inforf ${PLLMODULES_LIBRARIES} ortools)
target_include_directories(inforf PRIVATE ${LOCALDEPS_DIR}/include $<BUILD_INTERFACE:${ortools}>)

add_executable(rfdist main.cpp)
target_link_libraries(rfdist inforf)
target_include_directories(rfdist SYSTEM PUBLIC $<BUILD_INTERFACE:${LOCALDEPS_DIR}/include>)

add_subdirectory(io)
