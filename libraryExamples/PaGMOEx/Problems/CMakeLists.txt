 #    Copyright (c) 2010-2016, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(MY_PAGMO_PROBLEMS_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/earthMarsTransfer.cpp"
)

# Set the header files.
set(MY_PAGMO_PROBLEMS_HEADERS
  "${CMAKE_CURRENT_SOURCE_DIR}/earthMarsTransfer.h"
)

# Add static libraries, can be used later on outside this application.
add_library(my_pagmo_problems STATIC ${MY_PAGMO_PROBLEMS_SOURCES} ${MY_PAGMO_PROBLEMS_HEADERS})
