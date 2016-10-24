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
    using namespace propulsion;

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
    const double simulationEndEpoch = 5000.0 * tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 3600.0;

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Sun", "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects.
    double vehicleMass = 2000.0;
    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Asterix" ]->setConstantBodyMass( vehicleMass );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.

    boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionSettings =
            boost::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                "Sun", true, false );
    boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings =
            boost::make_shared< ConstantThrustEngineSettings >(
                5.0E-2, 4000.0 );

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity) );
    accelerationsOfAsterix[ "Asterix" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                       thrustDirectionSettings, thrustMagnitudeSettings ) );

    accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Sun" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Asterix" ] = createMassRateModel( "Asterix", boost::make_shared< FromThrustMassModelSettings >( 1 ),
                                                       bodyMap, accelerationModelMap );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Set initial conditions for satellites that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to
    // Cartesian elements.


    // Set initial state
    Eigen::VectorXd systemInitialState = spice_interface::getBodyCartesianStateAtEpoch(
                "Earth", "Sun", "ECLIPJ2000", "None", 0.0 );

    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch );

    boost::shared_ptr< PropagatorSettings< double > > massPropagatorSettings =
            boost::make_shared< MassPropagatorSettings< double > >(
                boost::assign::list_of( "Asterix" ), massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );


    std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Write Asterix propagation history to file.
    writeDataMapToTextFile( integrationResult,
                            "thrustPropagationHistory.dat",
                            tudat_applications::getOutputPath( ),
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;}
