 #    Copyright (c) 2010-2016, Delft University of Technology
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
project(Pagmo2LibraryExample)


# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../..")

if(NOT LIBROOT)
  set(LIBROOT "${CODEROOT}/lib")
endif()
if(NOT BINROOT)
  set(BINROOT "${CODEROOT}/bin")
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

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
             COMPONENTS thread date_time system unit_test_framework serialization filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Find PaGMO library on local system.
if( USE_PAGMO )
  find_package(PaGMO2)

  # Include PaGMO directories.
  if(NOT APPLE)
    include_directories(SYSTEM AFTER "${PAGMO_INCLUDE_DIR}")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${PAGMO_INCLUDE_DIR}\"")
  endif( )
endif( )

add_library(pagmo2_library_example STATIC "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/himmelblaus.hpp")

# Add application.
add_executable(application_PagmoHimmelblausExample "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

target_compile_options(application_PagmoHimmelblausExample PRIVATE "$<$<CONFIG:DEBUG>:${PAGMO_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${PAGMO_CXX_FLAGS_RELEASE}>")
set_property(TARGET application_PagmoHimmelblausExample PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
set_property(TARGET application_PagmoHimmelblausExample PROPERTY CXX_STANDARD 11)
set_property(TARGET application_PagmoHimmelblausExample PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET application_PagmoHimmelblausExample PROPERTY CXX_EXTENSIONS NO)
install(TARGETS application_PagmoHimmelblausExample RUNTIME DESTINATION "${BINROOT}/applications")

target_link_libraries(application_PagmoHimmelblausExample ${Boost_LIBRARIES} pthread)
