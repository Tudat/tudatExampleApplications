# - Try to find Tudat library
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Tudat 3.1.2)
# to require version 3.1.2 or newer of Tudat.
#
# Once done this will define
#
#  TUDAT_FOUND - system has tudat lib with correct version
#  TUDAT_INCLUDE_DIR - the tudat include directory
#
# This file is based on FindEigen3.cmake

# 
# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Copyright (c) 2012 Bryan Tong Minh <b.tongminh@student.tudelft.nl>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.

macro(_tudat_check_version)
  file(READ "${TUDAT_INCLUDE_DIR}/Tudat/tudatVersion.h" _tudat_header)

  string(REGEX MATCH "define[ \t]+TUDAT_VERSION_MAJOR[ \t]+([0-9]+)" _tudat_major_version_match "${_tudat_header}")
  set(TUDAT_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+TUDAT_VERSION_MINOR[ \t]+([0-9]+)" _tudat_minor_version_match "${_tudat_header}")
  set(TUDAT_MINOR_VERSION "${CMAKE_MATCH_1}")

  set(TUDAT_VERSION ${TUDAT_MAJOR_VERSION}.${TUDAT_MINOR_VERSION})
  if(${TUDAT_VERSION} VERSION_LESS ${Tudat_FIND_VERSION})
    set(TUDAT_VERSION_OK FALSE)
  else(${TUDAT_VERSION} VERSION_LESS ${Tudat_FIND_VERSION})
    set(TUDAT_VERSION_OK TRUE)
  endif(${TUDAT_VERSION} VERSION_LESS ${Tudat_FIND_VERSION})

  if(NOT TUDAT_VERSION_OK)

    message(STATUS "Tudat version ${TUDAT_VERSION} found in ${TUDAT_INCLUDE_DIR}, "
                   "but at least version ${Tudat_FIND_VERSION} is required")
  endif(NOT TUDAT_VERSION_OK)

  set(TUDAT_LIBRARIES "tudat")
  link_directories(${TUDAT_LIBRARIES_DIR})
endmacro(_tudat_check_version)

if (TUDAT_INCLUDE_DIR)

  # in cache already
  _tudat_check_version()
  set(TUDAT_FOUND ${TUDAT_VERSION_OK})

else (TUDAT_INCLUDE_DIR)

  find_path(TUDAT_BASE_PATH NAMES tudatVersion.h
      PATHS
      ${PROJECT_SOURCE_DIR}/External
      ${PROJECT_SOURCE_DIR}/../../../tudat/trunk
      ${PROJECT_SOURCE_DIR}/../../tudat/trunk
      ${PROJECT_SOURCE_DIR}/../../tudat
      ${CMAKE_INSTALL_PREFIX}/include
      PATH_SUFFIXES Tudat
    )
  set(TUDAT_INCLUDE_DIR ${TUDAT_BASE_PATH}/..)
  set(TUDAT_LIBRARIES_DIR ${TUDAT_BASE_PATH}/../lib)

  if(TUDAT_INCLUDE_DIR)
    _tudat_check_version()
  endif(TUDAT_INCLUDE_DIR)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Tudat DEFAULT_MSG TUDAT_INCLUDE_DIR TUDAT_VERSION_OK)

  mark_as_advanced(TUDAT_INCLUDE_DIR)

endif(TUDAT_INCLUDE_DIR)
