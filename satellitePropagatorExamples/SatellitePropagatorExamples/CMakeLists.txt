 #    Copyright (c) 2010-2013, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      120202    K. Kumar          Adapted to work with Earth orbiting satellite example.
 #      120523    A. Ronse          Updated terminology suiting application. Included script
 #                                  to find+load TudatCoreEnvironment.cmake.
 #      120205    K. Kumar          Removed obsolete TudatCoreEnvironment.cmake call; updated
 #                                  Tudat Core and Tudat minimum requirements; updated license;
 #                                  added static library so body.h shows up in project tree.
 #      130222    K. Kumar          Updated directory definitions; added unit test macro.
 #      130225    K. Kumar          Updated macros; updated set-function calls to check if
 #                                  variables have been set already; updated static library name
 #                                  and naming of applications.
 #     130326    K. Kumar           Added BUILD_DOCUMENTATION option for Doxygen html output
 #                                  generation; removed -O3 flag since it is set by "Release"
 #                                  build type if selected; updated install steps; re-organized
 #                                  file to isolate user-specified section.
 #
 #    References
 #
 #    Notes
 #

# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(SatellitePropagatorExamples)

# Load UserSettings.txt
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  message(STATUS "Building ${PROJECT_NAME} standalone.")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
else()
  message(STATUS "Building ${PROJECT_NAME} from within ${CMAKE_PROJECT_NAME}.")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  set(RELATIVE_PROJECT_PATH "${RELATIVE_PROJECT_PATH}" CACHE STRING "Relative path wrt to project for function")
  message(STATUS "Relative path (wrt to project): ${RELATIVE_PROJECT_PATH}")
endif()

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(STATUS "WARNING: building release version!")
    set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WARNING: building debug version!")
endif()

# Set module path to local CMake scripts.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(${CMAKE_MODULE_PATH}/compiler.cmake)

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
if(NOT APPLE)
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
find_package(Boost 1.45.0 
             COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

  # Find CSPICE library on local system.
  find_package(Spice)

  # Include CSpice directories.
  if(NOT APPLE)
    include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
  endif( )

# If documentation should be built, find Doxygen package and setup config file.
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)

  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it correctly.")
  endif()

  # Configure the Template Doxyfile for our specific project.
  configure_file(${CODEROOT}/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

  # Add a custom target to run Doxygen when ever the project is built.
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.
  add_custom_target (docs ALL COMMAND ${DOXYGEN_EXECUTABLE} 
                     ${PROJECT_BINARY_DIR}/Doxyfile 
                     SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)

  # Install documentation.
  install(DIRECTORY "${CODEROOT}/Documentation" 
          DESTINATION "${SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH}")
endif()


## Add single, Earth-orbiting satellite propagator application.
add_executable(application_SingleSatellitePropagator "${SRCROOT}/singleSatellitePropagator.cpp")
setup_executable_target(application_SingleSatellitePropagator "${SRCROOT}")
target_link_libraries(application_SingleSatellitePropagator
                      tudat_simulation_setup tudat_propagators tudat_aerodynamics tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES} )

add_executable(application_SinglePerturbedSatellitePropagator "${SRCROOT}/singlePerturbedSatellitePropagator.cpp")
setup_executable_target(application_SinglePerturbedSatellitePropagator "${SRCROOT}")
target_link_libraries(application_SinglePerturbedSatellitePropagator
                      tudat_simulation_setup tudat_propagators tudat_aerodynamics tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES} )

# Add Asterix and Obelix (satellites) propagator application.
add_executable(application_AsterixAndObelixPropagator "${SRCROOT}/asterixAndObelixPropagator.cpp")
setup_executable_target(application_AsterixAndObelixPropagator "${SRCROOT}")
target_link_libraries(application_AsterixAndObelixPropagator
                      tudat_simulation_setup tudat_propagators tudat_aerodynamics tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES} )

# Add Galileo constellation application.
add_executable(application_GalileoConstellationSimulator 
               "${SRCROOT}/galileoConstellationSimulator.cpp")
setup_executable_target(application_GalileoConstellationSimulator "${SRCROOT}")
target_link_libraries(application_GalileoConstellationSimulator
                      tudat_simulation_setup tudat_propagators tudat_aerodynamics tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES} )

# Add inner Solar system propagation
add_executable(application_InnerSolarSystemPropagation
               "${SRCROOT}/innerSolarSystemPropagation.cpp")
setup_executable_target(application_InnerSolarSystemPropagation "${SRCROOT}")
target_link_libraries(application_InnerSolarSystemPropagation
                      tudat_simulation_setup tudat_propagators tudat_aerodynamics tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES} )

# Add Apollo propagation
add_executable(application_ApolloEntry
               "${SRCROOT}/apolloCapsuleEntry.cpp")
setup_executable_target(application_ApolloEntry "${SRCROOT}")
target_link_libraries(application_ApolloEntry
                      tudat_simulation_setup tudat_propagators tudat_aerodynamics tudat_geometric_shapes tudat_gravitation tudat_mission_segments tudat_electro_magnetism tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_spice_interface cspice ${Boost_LIBRARIES} )
