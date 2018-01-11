# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(SatellitePropagatorExamples)

# Load UserSettings.txt
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(BUILD_STYLE "standalone")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
else()
  set(BUILD_STYLE "part of ${CMAKE_PROJECT_NAME}")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  set(RELATIVE_PROJECT_PATH "${RELATIVE_PROJECT_PATH}" CACHE STRING "Relative path wrt to project for function")
  # message(STATUS "Relative path (wrt to project): ${RELATIVE_PROJECT_PATH}")
endif()

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
  set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_BUILD_TYPE Debug)
endif()

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

# Add local module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")
message(STATUS "CMake Module path(s): ${CMAKE_MODULE_PATH}")

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(compiler)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Set testing options based on platform.
enable_testing()

# Set lib and bin directories where static libraries and unit tests are built.
if(NOT LIBROOT)
  set(LIBROOT "${CODEROOT}/lib")
endif()
if(NOT BINROOT)
  set(BINROOT "${CODEROOT}/bin")
endif()

# Set the global macros for setting up targets.
macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_target)

# Define the install targets to create a distribution.
if(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)
    set(TUDAT_BUNDLE_DISTRIBUTION_PATH "${CODEROOT}")
endif(NOT TUDAT_BUNDLE_DISTRIBUTION_PATH)

if(NOT SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH)
    set(SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH
        "${TUDAT_BUNDLE_DISTRIBUTION_PATH}/tudatApplications/satellitePropagatorExamples")
endif(NOT SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH)

# Install files.
install(DIRECTORY "${SRCROOT}/"
        DESTINATION "${SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH}/SatellitePropagatorExamples"
        PATTERN ".DS_STORE" EXCLUDE
        PATTERN "CMakeLists.txt.user" EXCLUDE
        PATTERN ".svn" EXCLUDE
        PATTERN ".git" EXCLUDE
        PATTERN ".bzr" EXCLUDE
)

# Include the top-level directories.
include_directories(AFTER
  "${CODEROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif()

# Configure Boost libraries.
if(NOT Boost_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
endif()
if(NOT Boost_USE_MULTITHREADED)
  set(Boost_USE_MULTITHREADED ON)
endif()
if(NOT Boost_USE_STATIC_RUNTIME)
  set(Boost_USE_STATIC_RUNTIME ON)
endif()

# Find Boost libraries on local system.
find_package(Boost 1.55.0
             COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE OR APPLE_INCLUDE_FORCE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

# Find PaGMO library on local system.
if( USE_OPTIMIZATION_PACKAGE )
  find_package(PaGMO2)

  # Include PaGMO directories.
  if(NOT APPLE)
    include_directories(SYSTEM AFTER "${PAGMO_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${PAGMO_INCLUDE_DIR}\"")
  endif( )

  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES pthread)

endif( )

find_package(NLOPT REQUIRED)
message(STATUS "NLopt include directory: ${NLOPT_INCLUDE_DIR}")
message(STATUS "NLopt library: ${NLOPT_LIBRARY}")
include_directories(SYSTEM AFTER "${NLOPT_INCLUDE_DIR}")

set( PAGMO_HIMMELBLAUS_FUNCTION_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/himmelblaus.hpp")

add_library(pagmo2_himmelblau_example_library "${PAGMO_HIMMELBLAUS_FUNCTION_HEADER}" "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
target_compile_options(pagmo2_himmelblau_example_library PRIVATE "$<$<CONFIG:DEBUG>:${PAGMO_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${PAGMO_CXX_FLAGS_RELEASE}>")
set_property(TARGET pagmo2_himmelblau_example_library PROPERTY CXX_STANDARD 11)
set_property(TARGET pagmo2_himmelblau_example_library PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET pagmo2_himmelblau_example_library PROPERTY CXX_EXTENSIONS NO)
install(TARGETS pagmo2_himmelblau_example_library DESTINATION "${BINROOT}/lib")

# Add application.
add_executable(application_PagmoHimmelblausExample "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

target_compile_options(application_PagmoHimmelblausExample PRIVATE "$<$<CONFIG:DEBUG>:${PAGMO_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${PAGMO_CXX_FLAGS_RELEASE}>")
set_property(TARGET application_PagmoHimmelblausExample PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
set_property(TARGET application_PagmoHimmelblausExample PROPERTY CXX_STANDARD 11)
set_property(TARGET application_PagmoHimmelblausExample PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET application_PagmoHimmelblausExample PROPERTY CXX_EXTENSIONS NO)
install(TARGETS application_PagmoHimmelblausExample RUNTIME DESTINATION "${BINROOT}/applications")

target_link_libraries(application_PagmoHimmelblausExample pthread ${Boost_LIBRARIES} )

add_executable(application_OptimizerComparison "${CMAKE_CURRENT_SOURCE_DIR}/optimizerComparison.cpp")

target_compile_options(application_OptimizerComparison PRIVATE "$<$<CONFIG:DEBUG>:${PAGMO_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${PAGMO_CXX_FLAGS_RELEASE}>")
set_property(TARGET application_OptimizerComparison PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
set_property(TARGET application_OptimizerComparison PROPERTY CXX_STANDARD 11)
set_property(TARGET application_OptimizerComparison PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET application_OptimizerComparison PROPERTY CXX_EXTENSIONS NO)
install(TARGETS application_OptimizerComparison RUNTIME DESTINATION "${BINROOT}/applications")

target_link_libraries(application_OptimizerComparison pthread nlopt ${Boost_LIBRARIES} )

add_executable(application_MultiObjectiveOptimizerComparison "${CMAKE_CURRENT_SOURCE_DIR}/multiObjectiveOptimizerComparison.cpp")

target_compile_options(application_MultiObjectiveOptimizerComparison PRIVATE "$<$<CONFIG:DEBUG>:${PAGMO_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${PAGMO_CXX_FLAGS_RELEASE}>")
set_property(TARGET application_MultiObjectiveOptimizerComparison PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
set_property(TARGET application_MultiObjectiveOptimizerComparison PROPERTY CXX_STANDARD 11)
set_property(TARGET application_MultiObjectiveOptimizerComparison PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET application_MultiObjectiveOptimizerComparison PROPERTY CXX_EXTENSIONS NO)
install(TARGETS application_MultiObjectiveOptimizerComparison RUNTIME DESTINATION "${BINROOT}/applications")

target_link_libraries(application_MultiObjectiveOptimizerComparison pthread nlopt ${Boost_LIBRARIES} )

