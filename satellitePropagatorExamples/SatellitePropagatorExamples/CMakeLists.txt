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
 #
 #    References
 #
 #    Notes
 #

# Specify minimum CMake version required to compile application.
cmake_minimum_required( VERSION 2.6 )

# Specific project name.
project( SatelitePropagatorExamples )

# Set path to local CMake modules.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the code root directory.
set(CODEROOT "${SRCROOT}/..")

# Set lib and bin directories where static libraries, applications, and unit tests are built.
set(LIBROOT "${CODEROOT}/lib")
set(BINROOT "${CODEROOT}/bin")

# Set testing options based on platform.
enable_testing()

# Set platform-specific compiler flags.
if (WIN32)
  if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
    # Because we are using static boost libraries, with static runtime, we need to force MSVC to
    # also use static runtime: (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace)
    foreach(flag_var
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      # Find all dynamic runtime (MD) references and replace with static (MT)
      if(${flag_var} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
endif()

# Defines the install target to create a distribution.
if(NOT SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH)
    set(SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH "${CODEROOT}/dist/satellitePropagatorExamples")
endif(NOT SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH)

# ------------------------------------------------------------------------------------
# Don't define any new variables with SET that the user might change beyond this point.
#
# Load user specific settings. These can overwrite all previous settings.
if(NOT SATELLITE_PROPAGATOR_EXAMPLES_USER_SETTINGS AND NOT DEFINED ENV{SATELLITE_PROPAGATOR_EXAMPLES_USER_SETTINGS})
  # Load from default paths, least important first.
  include( "../UserSettings.txt" OPTIONAL)
  include( "UserSettings.txt" OPTIONAL)
else()
  message(STATUS "Loading preset...")
  # Load the user specified file.
  if( SATELLITE_PROPAGATOR_EXAMPLES_USER_SETTINGS )
    include( $SATELLITE_PROPAGATOR_EXAMPLES_USER_SETTINGS OPTIONAL)
  else()
    include( $ENV{SATELLITE_PROPAGATOR_EXAMPLES_USER_SETTINGS} OPTIONAL)
  endif()
endif()

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release". NB: It is important this
# section is placed AFTER the "UserSettings.txt" file is included, to ensure that that the correct
# indiciation is given for which build type has been selected.
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(STATUS "WARNING: building release version!")
    set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WARNING: building debug version!")
endif()

# Install files.
install(DIRECTORY "${SRCROOT}/"
        DESTINATION "${SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH}/SatellitePropagatorExamples"
        FILES_MATCHING
            PATTERN "*"
            PATTERN "CMakeLists.txt.user" EXCLUDE
)

install(FILES ${CODEROOT}/LICENSE ${CODEROOT}/NOTICE ${CODEROOT}/README ${CODEROOT}/UserSettings.txt.example
  DESTINATION "${SATELLITE_PROPAGATOR_EXAMPLES_DISTRIBUTION_PATH}")

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
find_package(Boost 1.45.0 COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find Tudat Core library on local system.
find_package(TudatCore 2.0 REQUIRED)

# Include Tudat Core directories.
# Set CMake flag to suppress TudatCore warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_CORE_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_CORE_INCLUDE_DIR}\"")
endif( )

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

# Add library.
add_library(satellite_propagator_examples STATIC body.h void.cpp)
setup_library_target(satellite_propagator_examples "${SRCROOT}")

# Add Asterix and Obelix (satellites) propagator application.
add_executable(application_AsterixAndObelixPropagator "${SRCROOT}/asterixAndObelixPropagator.cpp")
setup_executable_target(application_AsterixAndObelixPropagator "${SRCROOT}")
target_link_libraries(application_AsterixAndObelixPropagator satellite_propagator_examples tudat_gravitation ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )

# Add Galileo constellation application.
add_executable(application_GalileoConstellationSimulator "${SRCROOT}/galileoConstellationSimulator.cpp")
setup_executable_target(application_GalileoConstellationSimulator "${SRCROOT}")
target_link_libraries(application_GalileoConstellationSimulator satellite_propagator_examples tudat_gravitation ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES} )