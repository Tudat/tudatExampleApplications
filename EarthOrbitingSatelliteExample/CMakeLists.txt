# \file CMakeLists.txt
 #    This is the main CMake file for the Earth-orbiting satellite example Tudat application.
 #
 #    Path              : /EarthOrbitingSatelliteExample/
 #    Version           : 2
 #    Check status      : Checked
 #
 #    Author/Checker    : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Author            : K. Kumar
 #    Affiliation       : Delft University of Technology
 #    E-mail address    : K.Kumar@tudelft.nl
 #
 #    Date created      : 20 August, 2011
 #    Last modified     : 2 February, 2012
 #
 #    References
 #
 #    Notes
 #
 #    Copyright (c) 2010-2011 Delft University of Technology.
 #
 #    This software is protected by national and international copyright.
 #    Any unauthorized use, reproduction or modification is unlawful and
 #    will be prosecuted. Commercial and non-private application of the
 #    software in any form is strictly prohibited unless otherwise granted
 #    by the authors.
 #
 #    The code is provided without any warranty; without even the implied
 #    warranty of merchantibility or fitness for a particular purpose.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #	120202    K. Kumar 	    Adapted to work with Tudat and Tudat Core.
 #/

# Specify minimum CMake version required to compile Tudat.
cmake_minimum_required( VERSION 2.6 )

# Specific Tudat project name.
project( EarthOrbitingSatelliteExample )

# Set module path to local CMake scripts.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)

# Set root-directory for code to current source directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the basic input (trunk) and output (lib & include) directories.
set(SRCROOT "${CODEROOT}")

# Set CMake build-type.
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)

# Set testing options based on platform.
enable_testing()

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
    set(LIBROOT "${SRCROOT}/lib")
    set(INCLUDEROOT "${SRCROOT}/include")
    set(BINROOT "${SRCROOT}/bin")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(LIBROOT "${SRCROOT}/lib")
    set(INCLUDEROOT "${SRCROOT}/include")
    set(BINROOT "${SRCROOT}/bin")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  set(LIBROOT "${SRCROOT}/lib")
  set(INCLUDEROOT "${SRCROOT}/include")
  set(BINROOT "${SRCROOT}/bin")
endif()


# Set the global macros for setting up targets.
macro(setup_tudat_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_tudat_executable_target)

macro(setup_tudat_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
  install(TARGETS ${target_name} LIBRARY DESTINATION "${BINROOT}" ARCHIVE DESTINATION "${LIBROOT}")
endmacro(setup_tudat_library_target)

macro(setup_custom_test_program target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_custom_test_program)

macro(setup_headers HEADER_FILES HEADER_PATH)
    foreach(CURRENT_HEADER_FILE ${HEADER_FILES})
      install(FILES "${CURRENT_HEADER_FILE}" DESTINATION "${INCLUDEROOT}${HEADER_PATH}")
    endforeach(CURRENT_HEADER_FILE)
endmacro(setup_headers)

# Include the top-level directories.
include_directories(AFTER
  "${SRCROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)
# If Eigen3 is not found, use the following to manually set the path to the Eigen3 directory.
# A default location is given as an example.
#set(EIGEN3_INCLUDE_DIR "${SRCROOT}/External/eigen3")

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif()

# Configure Boost libraries.
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)

# Find Boost libraries on local system.
find_package(Boost 1.45.0 COMPONENTS thread date_time system unit_test_framework filesystem regex iostreams REQUIRED)
# If Boost is not found, use the following to manually set the path to the Boost include and
# library directories.
#SET(Boost_INCLUDE_DIRS "C:\My\Custom\Path\To\Boost\Include\Directory")
#SET(Boost_LIBRARY_DIRS "C:\My\Custom\Path\To\Boost\Library\Directory")

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find Tudat Core library on local system.
find_package(TudatCore 0.1 REQUIRED)
# If TudatCore is not found, use the following to manually set the path to the TudatCore include and
# library directories.
#SET(TUDAT_CORE_INCLUDE_DIR "C:\My\Custom\Path\To\TudatCore\Include\Directory")
#SET(TUDAT_CORE_LIBRARIES_DIR "C:\My\Custom\Path\To\TudatCore\Library\Directory")

# Include TudatCore directories.
# Set CMake flag to suppress TudatCore warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_CORE_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_CORE_INCLUDE_DIR}\"")
endif()

# Find Tudat library on local system.
find_package(Tudat 0.1 REQUIRED)
# If Tudat is not found, use the following to manually set the path to the Tudat include and
# library directories.
#SET(TUDAT_INCLUDE_DIR "C:\My\Custom\Path\To\Tudat\Include\Directory")
#SET(TUDAT_LIBRARIES_DIR "C:\My\Custom\Path\To\Tudat\Library\Directory")

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

# Add Earth-orbiting satellite example application
add_executable(application_EarthOrbitingSatelliteExample "${SRCROOT}/earthOrbitingSatelliteExample.cpp")
setup_tudat_executable_target(application_EarthOrbitingSatelliteExample "${SRCROOT}")
target_link_libraries(application_EarthOrbitingSatelliteExample tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics tudat_core ${Boost_LIBRARIES})

# End of file.
