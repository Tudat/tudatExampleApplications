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
#include <fstream>

#include "Problems/earthMarsTransfer.h"


using namespace pagmo;

//! Execute  main
int main( )
{
   // Set the PRNG seed, such that results are reproducable
   int seed = 123;
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
   problem prob{EarthMarsTransfer( bounds )};

   // Select the self-adaptive differential evolution algorithm.
   // One generation per evolution step.
   algorithm algo{de1220(10)};

   // Create an island with 8 individuals
   island isl{algo, prob, 128};

   int i = 0;
   // For 30 generations optimise the population in the island
   for(  ; i < 30; i++ )
   {
        isl.evolve();
        while( isl.status()!=pagmo::evolve_status::idle )
            isl.wait();
        int c = isl.get_population().best_idx();
        vector_double cx = isl.get_population().champion_x();
        vector_double cf = isl.get_population().champion_f();
        print("GEN=", i, " ID=", c, " DV=", cf[ 0 ], "m/s DEP=", cx[ 0 ],
               "JD TOF=", cx[ 1 ], "d\n" );
   }

   return 0;

}
