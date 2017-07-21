/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>
/*#include <cmath>
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <Eigen/Core>
*/
#include <pagmo/pagmo.hpp>
#include "Pagmo2Example/Problems/earthMarsTransfer.h"

using namespace pagmo;


//! Execute main
int main( )
{
    // Set the PRNG seed, such that results are reproducable
    int seed = 123456;
    pagmo::random_device::set_seed( seed );

    // We have two decision variables each with a lower and upper
    // bound, create a vector of vectors that will contain these.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 2, 0.0 ) );

    // Search between 2020 and 2025 for flight duration between 200
    // and 1000 days.
    bounds[ 0 ][ 0 ] = 2458849.5;
    bounds[ 1 ][ 0 ] = 2460676.5;
    bounds[ 0 ][ 1 ] = 200;
    bounds[ 1 ][ 1 ] = 1000;

    // Define the problem
    problem prob(EarthMarsTransfer(bounds));

    // Select the self-adaptive differential evolution algorithm.
    // One generation per evolution step.
    algorithm algo{sade(1)};


    // Create an archipelago with 1 island with 8 individuals
    archipelago archi{1, algo, prob, 8};

    // For 25 generation optimise the population in the island
    for( int i=0 ; i < 25; i++ )
    {
        archi.evolve();
        int c = archi[0].get_population().best_idx();
        vector_double cx = archi[0].get_population().champion_x();
        vector_double  cf = archi[0].get_population().champion_f();
        std::cout << "GEN=" << i << " ID=" << c << " DV=" << cf[ 0 ]
                  << "m/s DEP=" << cx[ 0 ] << "JD TOF=" << cx[ 1 ] << "d" << std::endl;
    }

    return 0;
}
