/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Gravitation/centralJ2J3J4GravityModel.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include <Tudat/Astrodynamics/Propagators/dynamicsSimulator.h>
#include <Tudat/External/SpiceInterface/spiceInterface.h>
#include <Tudat/SimulationSetup/body.h>
#include <Tudat/SimulationSetup/createAccelerationModels.h>
#include <Tudat/SimulationSetup/defaultBodies.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <utility>

#include <Eigen/Core>

#include <SatellitePropagatorExamples/applicationOutput.h>

//! Execute propagation of orbits of Asterix and Obelix around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace basic_astrodynamics;
    using namespace basic_mathematics;
    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using namespace input_output;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLES      //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 60.0;

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );

    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
                basic_mathematics::Vector6d::Zero( ), "SSB", "J2000" );

    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = NULL;
    bodySettings[ "Earth" ]->shapeModelSettings = NULL;

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Obelix" ] = boost::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationMap[  "Asterix" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfObelix;
    accelerationsOfObelix[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationMap[  "Obelix" ] = accelerationsOfObelix;
    bodiesToPropagate.push_back( "Obelix" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Set initial conditions for satellites that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to
    // Cartesian elements.

    // Set Keplerian elements for Asterix.
    Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0e3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Set Keplerian elements for Obelix.
    Vector6d obelixInitialStateInKeplerianElements( 6 );
    obelixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 12040.6e3;
    obelixInitialStateInKeplerianElements( eccentricityIndex ) = 0.4;
    obelixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( -23.5 );
    obelixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 10.6 );
    obelixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 367.9 );
    obelixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 93.4 );

    // Convert initial states from Keplerian to Cartesian elements.

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    const Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );

    // Convert Obelix state from Keplerian elements to Cartesian elements.
    const Vector6d obelixInitialState = convertKeplerianToCartesianElements(
                obelixInitialStateInKeplerianElements,
                earthGravitationalParameter );


    // Set initial state
    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 12 );
    systemInitialState.segment( 0, 6 ) = asterixInitialState;
    systemInitialState.segment( 6, 6 ) = obelixInitialState;

    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, simulationEndEpoch, fixedStepSize );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Retrieve numerically integrated states of vehicles.
    std::map< double, Eigen::VectorXd > asterixPropagationHistory;
    std::map< double, Eigen::VectorXd > obelixPropagationHistory;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
         stateIterator != integrationResult.end( ); stateIterator++ )
    {
        asterixPropagationHistory[ stateIterator->first ] = stateIterator->second.segment( 0, 6 );
        obelixPropagationHistory[ stateIterator->first ] = stateIterator->second.segment( 6, 6 );
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Write Asterix propagation history to file.
    writeDataMapToTextFile( asterixPropagationHistory,
                            "asterixPropagationHistory.dat",
                            tudat_applications::getOutputPath( ),
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    // Write obelix propagation history to file.
    writeDataMapToTextFile( obelixPropagationHistory,
                            "obelixPropagationHistory.dat",
                            tudat_applications::getOutputPath( ),
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;}
