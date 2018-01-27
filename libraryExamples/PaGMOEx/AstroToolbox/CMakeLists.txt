#     Copyright (c) 2010-2012 Delft University of Technology.
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
#      111025    K. Kumar          Adapted file to work with Revision 194.
#      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
#      130111    D. Dirkx          Modified for separate spice compilation
#      151111    J. Geul           Make file more targeted for CSPICE
#/

# Specify minimum CMake version required to compile cspice.
cmake_minimum_required( VERSION 2.6 )

# Specific project name.
project( AstroToolbox )

# Detect if CSPICE is build standalone or from within project
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  message(STATUS "Building ${PROJECT_NAME} standalone.")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
else()
  message(STATUS "Building ${PROJECT_NAME} from within ${CMAKE_PROJECT_NAME}.")
  include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
  set(PROJECTROOT ${PROJECTROOT} PARENT_SCOPE)
endif()

# Define the basic input (trunk).
set(SRCROOT ${PROJECT_SOURCE_DIR})

set(LIBROOT "${SRCROOT}/lib")
set(BINROOT "${SRCROOT}/bin")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -Wno-incompatible-pointer-types")
# Hide warnings for implicit-int and incompatible-pointer-types

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

# Set the global macros for setting up targets.
macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
 set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
 set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
 install(TARGETS ${target_name} LIBRARY DESTINATION "${BINROOT}" ARCHIVE DESTINATION "${LIBROOT}")
endmacro(setup_library_target)


file(GLOB_RECURSE ASTRO_TOOLBOX_HEADERS ${SRCROOT} RELATIVE ${SRCROOT} *.h)
file(GLOB_RECURSE ASTRO_TOOLBOX_SOURCES ${SRCROOT} RELATIVE ${SRCROOT} *.cpp)


# Add static libraries.
add_library(astro_toolbox STATIC ${ASTRO_TOOLBOX_SOURCES} ${ASTRO_TOOLBOX_HEADERS})
setup_library_target(astro_toolbox "${SRCROOT}")

# End of file.
