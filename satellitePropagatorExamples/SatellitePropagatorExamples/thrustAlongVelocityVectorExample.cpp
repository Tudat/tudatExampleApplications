
/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "SatellitePropagatorExamples/applicationOutput.h"


//! Execute propagation of orbit of Asterix around the Earth.
int main()
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_mathematics;
    using namespace tudat::basic_astrodynamics;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Load spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Create Earth object
    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );

    NamedBodyMap bodyMap = createBodies( bodySettings );
    // Create vehicle objects.
    double vehicleMass = 5.0E3;
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define thrust settings
    double thrustMagnitude = 25.0;
    double specificImpulse = 5000.0;
    boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            boost::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                "Earth", true, false );
    boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings =
            boost::make_shared< ConstantThrustEngineSettings >(
                thrustMagnitude, specificImpulse );


    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Vehicle" ].push_back(
                boost::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings) );
    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );


    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set initial state
    Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );
    systemInitialState( 0 ) = 8.0E6;
    systemInitialState( 4 ) = 7.5E3;

    // Define propagation termination conditions (stop after 2 weeks).
    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( 14.0 * physical_constants::JULIAN_DAY );

    // Define settings for propagation of translational dynamics.
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings,
              cowell );

    // Crete mass rate models
    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = createMassRateModel( "Vehicle", boost::make_shared< FromThrustMassModelSettings >( 1 ),
                                                      bodyMap, accelerationModelMap );

    // Create settings for propagating the mass of the vehicle
    boost::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            boost::make_shared< MassPropagatorSettings< double > >(
                boost::assign::list_of( "Vehicle" ), massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
                terminationSettings );

    // Create list of propagation settings.
    std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    // Create propagation settings for mass and translational dynamics concurrently
    boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, terminationSettings );

    // Define integrator settings
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 30.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );

    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Write satellite propagation history to file.
    input_output::writeDataMapToTextFile( numericalSolution,
                                          "velocityVectorThrustExample.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

}
