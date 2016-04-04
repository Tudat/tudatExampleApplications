/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>
#include <cmath>
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <Eigen/Core>

#include <pagmo/src/pagmo.h>
#include <pagmo/src/rng.h>

#include "PaGMOEx/Problems/earthMarsTransfer.h"

using namespace pagmo;
using boost::format;

//! Execute main
int main( )
{
    // Set the PRNG seed, such that results are reproducable
    int seed = 123456;
    pagmo::rng_generator::set_seed( seed );

    // If we have archipelagos, we also set the seed there
    // arch.set_seeds(sim_id);
    // Similarly set the seed for any other PRNG we might be using:
    // srand( seed );

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
    problem::EarthMarsTransfer prob( bounds );

    // Create a population (8 is minimum for jDE)
    population pop( prob, 8 );

    // Select the self-adaptive differential evolution algorithm.
    // One generation per evolution step.
    algorithm::jde algo( 1 );

    unsigned int i = 0;

    // For 25 generation optimise the population
    for( ; i < 25; i++ )
    {
        algo.evolve( pop );
        int c = pop.get_best_idx( );
        decision_vector cx = pop.get_individual( c ).cur_x;
        fitness_vector  cf = pop.get_individual( c ).cur_f;
        std::cout << "GEN=" << i << " ID=" << c << " DV=" << cf[ 0 ]
                  << "m/s DEP=" << cx[ 0 ] << "JD TOF=" << cx[ 1 ] << "d" << std::endl;
    }

    // Demonstrate serialization
    // Serialization is the process of converting objects (like
    // instances classes) into a stream of bytes. This allows us to
    // save a population disk and reload it at another time (to
    // continue the optimization).
    std::string appPath( __FILE__ );
    appPath = appPath.substr( 0, appPath.find_last_of( "/\\") + 1 );

    // Save the population (contains also the problem) and some other
    // variable (just to show how its done).
    std::cout << "Saving the population" << std::endl;
    std::ofstream ofs( appPath + "pop.bak" );
    boost::archive::text_oarchive oa( ofs );
    oa << i << bounds << pop;
    ofs.close();

    // Change the value of "i" and the population
    i = 999;
    pop.clear( );
    std::cout << "Clearing i=" << i << " Npop=" << pop.size( ) << std::endl;

    // Reload everything ("i" will reset to the stored value)
    std::cout << "Loading the population" << std::endl;
    std::ifstream ifs( appPath + "pop.bak" );
    boost::archive::text_iarchive ia( ifs );
    ia >> i >> bounds >> pop;
    ifs.close( );

    // Resume the evolution
    for( ; i < 70; i++ )
    {
        algo.evolve( pop );
        int c = pop.get_best_idx( );
        decision_vector cx = pop.get_individual( c ).cur_x;
        fitness_vector  cf = pop.get_individual( c ).cur_f;
        std::cout << "GEN=" << i << " ID=" << c << " DV=" << cf[ 0 ]
                  << "m/s DEP=" << cx[ 0 ] << "JD TOF=" << cx[ 1 ] << "d" << std::endl;
    }

    std::cout << "How great that we were able to load the population and continue"
              << " where we left off!" << std::endl;

    return 0;
}
