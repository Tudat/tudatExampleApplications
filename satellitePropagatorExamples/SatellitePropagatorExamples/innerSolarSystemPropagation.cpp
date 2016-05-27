/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Propagators/nBodyCowellStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsSimulator.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/SimulationSetup/createBodies.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/defaultBodies.h"

#include "SatellitePropagatorExamples/applicationOutput.h"

//Using declarations.
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat;

int main( )
{

    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");

    // Define bodies in simulation.
    unsigned int totalNumberOfBodies = 6;
    std::vector< std::string > bodyNames;
    bodyNames.resize( totalNumberOfBodies );
    bodyNames[ 0 ] = "Moon";
    bodyNames[ 1 ] = "Earth";
    bodyNames[ 2 ] = "Mars";
    bodyNames[ 3 ] = "Venus";
    bodyNames[ 4 ] = "Mercury";
    bodyNames[ 5 ] = "Sun";



    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 1.0E7 + 5.0 * physical_constants::JULIAN_YEAR;
    double maximumTimeStep = 3600.0;
    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    for( int centralBodySettings = 0; centralBodySettings < 2; centralBodySettings++ )
    {
        // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
        SelectedAccelerationMap accelerationMap;
        for( unsigned int i = 0; i < bodyNames.size( ); i++ )
        {
            std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > currentAccelerations;
            for( unsigned int j = 0; j < bodyNames.size( ); j++ )
            {
                if( i != j )
                {
                    currentAccelerations[ bodyNames.at( j ) ].push_back(
                                boost::make_shared< AccelerationSettings >( central_gravity ) );\
                }
            }
            accelerationMap[ bodyNames.at( i ) ] = currentAccelerations;
        }

        // Define list of bodies to propagate
        std::vector< std::string > bodiesToPropagate = bodyNames;
        unsigned int numberOfNumericalBodies = bodiesToPropagate.size( );

        // Define numerical integrator settings.
        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, initialEphemerisTime, finalEphemerisTime, 3600.0 );

        // Define central bodies to use in propagation.
        std::vector< std::string > centralBodies;
        centralBodies.resize( numberOfNumericalBodies );
        if( centralBodySettings == 0 )
        {
            for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
            {
                centralBodies[ i ] = "SSB";
            }
        }
        else if( centralBodySettings == 1 )
        {
            for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
            {
                if( i == 0 )
                {
                    centralBodies[ i ] = "Earth";
                }
                else if( i == 5 )
                {
                    centralBodies[ i ] = "SSB";
                }
                else
                {
                    centralBodies[ i ] = "Sun";
                }
            }
        }


        // Get initial state vector as input to integration.
        Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                    bodiesToPropagate, centralBodies, bodyMap, initialEphemerisTime );

        // Create acceleration models and propagation settings.
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );
        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );

        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        // Retrieve numerically integrated state for each body.
        std::vector< std::map< double, Eigen::VectorXd > > allBodiesPropagationHistory;
        allBodiesPropagationHistory.resize( numberOfNumericalBodies );
        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
             stateIterator != integrationResult.end( ); stateIterator++ )
        {
            for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
            {
                allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
            }
        }

        for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
        {
            // Write propagation history to file.
            input_output::writeDataMapToTextFile(
                        dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                        "innerSolarSystemPropagationHistory" + bodyNames.at( i ) +
                        boost::lexical_cast< std::string >( centralBodySettings ) + ".dat",
                        tudat_applications::getOutputPath( ),
                        "",
                        std::numeric_limits< double >::digits10,
                        std::numeric_limits< double >::digits10,
                        "," );
        }
    }
}

