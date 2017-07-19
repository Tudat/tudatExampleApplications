/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H
#define TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H

#include <vector>
#include <utility>

#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

#include <pagmo/pagmo.hpp>

#include <Eigen/Core>

typedef Eigen::Matrix< double, 6, 1 > StateType;

using namespace pagmo;

//! Test function for a new interplanetary trajectory class in Tudat
struct EarthMarsTransfer
{

    EarthMarsTransfer( const std::vector< std::vector< double > > problemBounds );

    // Calculates the fitness
    vector_double fitness( const vector_double &x ) const;

    std::pair<vector_double, vector_double> get_bounds() const;

    std::string get_name( ) const;

private:

    const std::vector< std::vector< double > > problemBounds_;

    StateType getPlanetPosition( const double date, const std::string planetName ) const;
};

//PAGMO_REGISTER_PROBLEM(EarthMarsTransfer)

#endif // TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H
