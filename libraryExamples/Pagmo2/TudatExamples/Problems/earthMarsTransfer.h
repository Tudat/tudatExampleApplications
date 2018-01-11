/*    Copyright (c) 2010-2017, Delft University of Technology
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
#include <limits>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include <Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h>

#include "pagmo/island.hpp"
//#include "pagmo/external/cereal/cereal.hpp"
#include "pagmo/io.hpp"
#include "pagmo/serialization.hpp"
#include "pagmo/problem.hpp"

#include <Eigen/Core>

typedef Eigen::Matrix< double, 6, 1 > StateType;

using namespace pagmo;

//! Test function for a new interplanetary trajectory class in Tudat
struct EarthMarsTransfer
{

    EarthMarsTransfer( const bool useTripTime = false ): useTripTime_( useTripTime ){ }

    EarthMarsTransfer( std::vector< std::vector< double > > &bounds, const bool useTripTime = false );

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    std::string get_name( ) const;

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        if(useTripTime_ )
        {
            return 2u;
        }
        else
        {
            return 1u;
        }

    }

private:

    const std::vector< std::vector< double > > problemBounds_;

    StateType getPlanetPosition( const double date, const std::string planetName ) const;

    bool useTripTime_;
};

#endif // TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H
