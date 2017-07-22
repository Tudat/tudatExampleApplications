 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the header files.
set(MY_PAGMO_PROBLEMS_HEADERS
  "${CMAKE_CURRENT_SOURCE_DIR}/easy_problem.hpp"
)

list(APPEND SOURCE_FILES ${MY_PAGMO_PROBLEMS_HEADERS})
