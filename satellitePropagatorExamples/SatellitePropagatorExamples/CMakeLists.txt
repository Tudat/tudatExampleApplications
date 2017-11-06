 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

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
find_package(Boost 1.45.0 
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

  # Find CSPICE library on local system.
  find_package(Spice)

  # Include CSpice directories.
  if(NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
  endif( )

option(USE_NRLMSISE00 "build Tudat with NRLMSISE-00 enabled" ON)
if(NOT USE_NRLMSISE00)
  message(STATUS "NRLMSISE-00 disabled!")
  add_definitions(-DUSE_NRLMSISE00=0)
else()
  message(STATUS "NRLMSISE-00 enabled!")
  add_definitions(-DUSE_NRLMSISE00=1)
  # Find USE_NRLMSISE00 library on local system.
  find_package(NRLMSISE00)

  # Include NRLMSISE00 directories.
  if(NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${NRLMSISE00_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${NRLMSISE00_INCLUDE_DIR}\"")
  endif( )
endif( )

list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES "")
list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES "")
list(APPEND TUDAT_APPLICATION_ITRS_LIBRARIES "")


if(USE_SOFA)
 list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES sofa)
 list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES tudat_sofa_interface )
 list(APPEND TUDAT_APPLICATION_ITRS_LIBRARIES tudat_earth_orientation )
endif()

if(USE_CSPICE)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES cspice)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES tudat_spice_interface )
endif()

if(USE_JSON)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES nlohmann_json)
endif()

if(USE_NRLMSISE00)
  list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES nrlmsise00)
endif()


list(APPEND TUDAT_APPLICATION_PROPAGATION_LIBRARIES tudat_simulation_setup tudat_ground_stations tudat_propagators
    tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
    tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_APPLICATION_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
     tudat_statistics tudat_propagators ${TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_basic_mathematics
     tudat_input_output tudat_basics ${TUDAT_APPLICATION_EXTERNAL_LIBRARIES})

list(APPEND TUDAT_APPLICATION_ESTIMATION_LIBRARIES tudat_simulation_setup tudat_observation_models tudat_ground_stations tudat_acceleration_partials
    tudat_observation_partials tudat_estimatable_parameters tudat_orbit_determination  tudat_propagators
    tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
    tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_APPLICATION_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
    tudat_statistics tudat_propagators ${TUDAT_APPLICATION_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_basic_mathematics
    tudat_input_output tudat_basics ${TUDAT_APPLICATION_EXTERNAL_LIBRARIES})



## Add single, Earth-orbiting satellite propagator application.
add_executable(application_SingleSatellitePropagator "${SRCROOT}/singleSatellitePropagator.cpp")
setup_executable_target(application_SingleSatellitePropagator "${SRCROOT}")
target_link_libraries(application_SingleSatellitePropagator ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(application_SinglePerturbedSatellitePropagator "${SRCROOT}/singlePerturbedSatellitePropagator.cpp")
setup_executable_target(application_SinglePerturbedSatellitePropagator "${SRCROOT}")
target_link_libraries(application_SinglePerturbedSatellitePropagator
		      ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(application_VariationalEquationsPropagator "${SRCROOT}/variationalEquationsPropagatorExample.cpp")
setup_executable_target(application_VariationalEquationsPropagator "${SRCROOT}")
target_link_libraries(application_VariationalEquationsPropagator
                      ${TUDAT_APPLICATION_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add Asterix and Obelix (satellites) propagator application.
add_executable(application_AsterixAndObelixPropagator "${SRCROOT}/asterixAndObelixPropagator.cpp")
setup_executable_target(application_AsterixAndObelixPropagator "${SRCROOT}")
target_link_libraries(application_AsterixAndObelixPropagator
                      ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add Galileo constellation application.
add_executable(application_GalileoConstellationSimulator 
               "${SRCROOT}/galileoConstellationSimulator.cpp")
setup_executable_target(application_GalileoConstellationSimulator "${SRCROOT}")
target_link_libraries(application_GalileoConstellationSimulator
                      ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add inner Solar system propagation
add_executable(application_InnerSolarSystemPropagation
               "${SRCROOT}/innerSolarSystemPropagation.cpp")
setup_executable_target(application_InnerSolarSystemPropagation "${SRCROOT}")
target_link_libraries(application_InnerSolarSystemPropagation
                      ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add Apollo propagation
add_executable(application_ApolloEntry
               "${SRCROOT}/apolloCapsuleEntry.cpp")
setup_executable_target(application_ApolloEntry "${SRCROOT}")
target_link_libraries(application_ApolloEntry ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add JSON-based Apollo propagation
add_executable(application_ApolloEntryJSON
	       "${SRCROOT}/apolloCapsuleEntryJSON.cpp")
setup_executable_target(application_ApolloEntryJSON "${SRCROOT}")
target_link_libraries(application_ApolloEntryJSON json_interface_library ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add lifetime maximisation propagation
add_executable(application_LifetimeMaximisation
	       "${SRCROOT}/lifetimeMaximisation.cpp")
setup_executable_target(application_LifetimeMaximisation "${SRCROOT}")
target_link_libraries(application_LifetimeMaximisation json_interface_library ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add Thrust example 1
add_executable(application_ThrustAlongVelocityVector
               "${SRCROOT}/thrustAlongVelocityVectorExample.cpp")
setup_executable_target(application_ThrustAlongVelocityVector "${SRCROOT}")
target_link_libraries(application_ThrustAlongVelocityVector ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add Thrust example 2
add_executable(application_ThrustFromFile
               "${SRCROOT}/thrustAccelerationFromFileExample.cpp")
setup_executable_target(application_ThrustFromFile "${SRCROOT}")
target_link_libraries(application_ThrustFromFile ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# Add simulated Earth orbiter simulated POD example
add_executable(application_EarthOrbiterStateEstimation
               "${SRCROOT}/earthOrbiterStateEstimation.cpp")
setup_executable_target(application_EarthOrbiterStateEstimation "${SRCROOT}")
target_link_libraries(application_EarthOrbiterStateEstimation ${TUDAT_APPLICATION_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(application_EarthOrbiterBasicStateEstimation
               "${SRCROOT}/earthOrbiterBasicStateEstimation.cpp")
setup_executable_target(application_EarthOrbiterBasicStateEstimation "${SRCROOT}")
target_link_libraries(application_EarthOrbiterBasicStateEstimation ${TUDAT_APPLICATION_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )
