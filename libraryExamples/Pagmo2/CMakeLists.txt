project(Pagmo2Example)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(compiler)

find_package(Boost 1.55.0
    COMPONENTS thread date_time system serialization)
include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")

find_package(PaGMO)
include_directories(SYSTEM AFTER ${PAGMO_INCLUDE_DIR})

add_executable(pagmo_simple_problem pagmo_simple_problem.cpp)
target_link_libraries(pagmo_simple_problem pthread)
