add_library(inforf STATIC
    datastructures/PllTree.cpp
    datastructures/PllTree.hpp
    datastructures/PllSplits.cpp
    datastructures/PllSplits.hpp
    datastructures/SimilarityCache.cpp
    datastructures/SimilarityCache.hpp
    datastructures/UniquePllMap.cpp
    datastructures/UniquePllMap.hpp
    PhylogeneticMathUtils.cpp
    PhylogeneticMathUtils.hpp
    Metric.cpp
    Metric.hpp

)


set_target_properties(inforf PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

target_link_libraries(inforf Threads::Threads)
target_link_libraries(inforf ${PLLMODULES_LIBRARIES} ortools json_lib)
target_include_directories(inforf PRIVATE ${LOCALDEPS_DIR}/include $<BUILD_INTERFACE:${ortools}>)



add_library(json_lib INTERFACE)
message(${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(json_lib SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/io/json>)



add_executable(rfdist main.cpp)
#target_include_directories(rfdist PUBLIC rf_io_lib)
#target_link_libraries(rfdist rf_io_lib)
#target_link_libraries(rfdist json_lib)
target_link_libraries(rfdist inforf)
#target_link_libraries(rfdist gmp :libgmp.a)
target_include_directories(rfdist SYSTEM PUBLIC $<BUILD_INTERFACE:${LOCALDEPS_DIR}/include>)


add_subdirectory(io)
add_subdirectory(datastructures)
