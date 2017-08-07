/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/External/JsonInterface/simulation.h>

//! Execute propagation of orbit of Asterix around the Earth.
int main( int argc, char* argv[] )
{
    if ( argc != 2 )
    {
        std::cerr << "Usage: tudat \"relative or absolute path to a JSON input file\"" << std::endl;
        return EXIT_FAILURE;
    }
    else
    {
        std::string inputFilePath = argv[ 1 ];

        tudat::json_interface::Simulation< > simulation( inputFilePath );

        // simulation.integratorSettings->integratorType_ = tudat::numerical_integrators::rungeKuttaVariableStepSize;

        simulation.exportAsJSON( );

        simulation.run( );
        simulation.exportResults( );

        return EXIT_SUCCESS;
    }
}
