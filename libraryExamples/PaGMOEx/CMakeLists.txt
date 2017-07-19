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
project(Pagmo2LibraryExamples)

macro(add_tudat_application target_name main_path)
    add_executable("${target_name}" "${main_path}")
    target_link_libraries("${target_name}" pthread ${Boost_LIBRARIES})
endmacro(add_tudat_application)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")

include("${CMAKE_MODULE_PATH}/compiler.cmake")

set(LIBROOT "${CMAKE_CURRENT_SOURCE_DIR}/../lib")

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")

find_package(Boost 1.55.0
    COMPONENTS thread date_time system serialization)
include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")

find_package(PaGMO)
include_directories(SYSTEM AFTER ${PAGMO_INCLUDE_DIR})

find_package(Tudat 2.0 REQUIRED)
include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
include_directories(SYSTEM AFTER "${CMAKE_CURRENT_SOURCE_DIR}/..")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Problems")

add_tudat_application(pagmo_example_application "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
target_link_libraries(pagmo_example_application my_pagmo_problems
    tudat_basic_astrodynamics tudat_mission_segments)
