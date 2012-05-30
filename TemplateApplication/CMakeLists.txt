 #    Copyright (c) 2010-2012 Delft University of Technology.
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
 #	120202    K. Kumar 	    Adapted to work with Earth orbiting satellite example.
 #	120427	  A. Ronse	    Added option for using UserSettings.txt.
 #	120522	  A. Ronse	    Adapted for automatically finding TudatCoreEnvironment.cmake
 #			            	at standard location relative to project source.
 #
 #    References
 #
 #/

# Specify minimum CMake version required to compile application.
cmake_minimum_required( VERSION 2.6 )

# Specify application project name.
project( TemplateApplication )

# Set path to local CMake modules.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# Define the basic input (trunk) and output (lib & include) directories.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Set CMake build-type.
set(CMAKE_BUILD_TYPE Release)

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

# Set directories for compiled libraries and binaries.
set(LIBROOT "${SRCROOT}/lib")
set(BINROOT "${SRCROOT}/bin")

# Set the global macros for setting up targets. Do not modify!
# Different sub-directories are created for applications and unit tests.
macro(setup_application_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_application_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
  install(TARGETS ${target_name} LIBRARY DESTINATION "${BINROOT}" ARCHIVE DESTINATION "${LIBROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_executable_target)

# Include the top-level directories.
include_directories(AFTER
  "${SRCROOT}"
)

# Look for TudatCoreEnvironment.cmake in standard location.
if (NOT TUDAT_CORE_ENVIRONMENT_PATH)
  find_path(TUDAT_CORE_ENVIRONMENT_PATH TudatCoreEnvironment.cmake
      ${PROJECT_SOURCE_DIR}/../../tudatcore
      ${PROJECT_SOURCE_DIR}/../../../tudatcore/trunk
    )
endif()

# Load the environment file from the core distribution if specified or found. This sets the paths for TudatCore, Boost and Eigen.
if(TUDAT_CORE_ENVIRONMENT_PATH)
    message(STATUS "Loading Tudat Core distribution from ${TUDAT_CORE_ENVIRONMENT_PATH}/TudatCoreEnvironment.cmake")
    include("${TUDAT_CORE_ENVIRONMENT_PATH}/TudatCoreEnvironment.cmake")
endif(TUDAT_CORE_ENVIRONMENT_PATH)

# Load user specific settings. These can overwrite all previous settings. 
if(NOT TUDAT_USER_SETTINGS AND NOT DEFINED ENV{TUDAT_USER_SETTINGS})
  # Load from default paths, least important first.
  INCLUDE( "../UserSettings.txt" OPTIONAL)
  INCLUDE( "UserSettings.txt" OPTIONAL)
else()
	MESSAGE(STATUS "Loading custom preset user settings")
  # Load the user specified file.
  if( TUDAT_USER_SETTINGS )
    #INCLUDE( $TUDAT_USER_SETTINGS OPTIONAL)
  else()  
    #INCLUDE( $ENV{TUDAT_USER_SETTINGS} OPTIONAL)
  endif()
endif()

# Find Tudat Core library on local system.
find_package(TudatCore 1.1 REQUIRED)

# Include TudatCore directories.
# Set CMake flag to suppress TudatCore warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_CORE_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_CORE_INCLUDE_DIR}\"")
endif( )

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

# Find Tudat library on local system.
find_package(Tudat 1.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif()

# Add helloWorld application.
add_executable(application_HelloWorld "${SRCROOT}/helloWorld")
setup_application_executable_target(application_HelloWorld "${SRCROOT}")
target_link_libraries(application_HelloWorld tudat_gravitation ${TUDAT_CORE_LIBRARIES} ${Eigen_LIBRARIES} ${Boost_LIBRARIES} )
