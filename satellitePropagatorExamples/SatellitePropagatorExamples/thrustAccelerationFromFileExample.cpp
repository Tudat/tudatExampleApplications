/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/applicationOutput.h>


std::map< double, Eigen::Vector3d > getThrustData( )
{
    // Find filepath and folder of this cpp file
    std::string cppFilePath( __FILE__ );
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );

    // Load data into matrix
    Eigen::MatrixXd thrustForceMatrix =
            tudat::input_output::readMatrixFromFile( cppFolder + "testThrustValues.txt" , " \t" );

    // Fill thrustData map using thrustForceMatrix Eigen matrix
    std::map< double, Eigen::Vector3d > thrustData;
    for ( int i = 0; i < thrustForceMatrix.rows( ); i++ )
    {
        thrustData[ thrustForceMatrix( i, 0 ) ] = thrustForceMatrix.block( i, 1, 1, 3 ).transpose( );
    }
    return thrustData;
}


//! Execute propagation of orbit of vehicle around the Earth. The vehicle is subject to a thrustforce, which is specified in
//! the nonconstantThrust.txt file. In that file, the first column is time in seconds, the last three columns give the x, y
//! and z components of the thrust force in the J2000 (?) frame.
int main()
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::numerical_integrators;
    using namespace tudat::interpolators;
    using namespace tudat::unit_conversions;
    using namespace tudat::propagators;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Define body settings for simulation.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    double vehicleMass = 2000.0;
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Retrieve thrust data as function of time.
    std::map< double, Eigen::Vector3d > thrustData = getThrustData( );

    // Make interpolator
    boost::shared_ptr< InterpolatorSettings >
            thrustInterpolatorSettingsPointer = boost::make_shared< InterpolatorSettings >( linear_interpolator );

    // Creating settings for thrust force
    boost::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector3d > >
            thrustInterpolatorPointer = createOneDimensionalInterpolator< double, Eigen::Vector3d >(
                thrustData, thrustInterpolatorSettingsPointer );

    // Define specific impulse
    double constantSpecificImpulse = 3000.0;

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Vehicle" ].push_back(
                boost::make_shared< ThrustAccelerationSettings >(
                    thrustInterpolatorPointer,
                    boost::lambda::constant( constantSpecificImpulse ), lvlh_thrust_frame, "Earth" ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Set initial conditions for the vehicle satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.
    // Set Keplerian elements for vehicle.
    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 72130.0e3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.6;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 169.0 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 45.0 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 80.0 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 15.0 );

    // Convert vehicle state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                vehicleInitialStateInKeplerianElements,
                earthGravitationalParameter );

    // Define propagation termination conditions (stop after 2 weeks).
    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( 4.0E6 );

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

    // Define list of dependent variables to save.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    lvlh_to_inertial_frame_rotation_dependent_variable, "Vehicle", "Earth" ) );

    // Create object with list of dependent variables
    boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    // Create propagation settings for mass and translational dynamics concurrently
    boost::shared_ptr< PropagatorSettings< > > propagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

    // Define integrator settings
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 30.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT AND PRINT OUTPUT TO FILE         //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

    // Manually add thrust force in LVLH frame to output
    for( std::map< double, Eigen::VectorXd >::iterator outputIterator = dependentVariableResult.begin( );
         outputIterator != dependentVariableResult.end( ); outputIterator++ )
    {
        Eigen::Matrix3d currentRotationMatrix =
                getMatrixFromVectorRotationRepresentation( outputIterator->second.segment( 3, 9 ) );
        Eigen::Vector3d currentThrust = outputIterator->second.segment( 0, 3 );
        Eigen::VectorXd newOutput = Eigen::VectorXd( 15 );
        newOutput.segment( 0, 12 ) = outputIterator->second;
        newOutput.segment( 12, 3 ) =
                integrationResult.at( outputIterator->first )( 6 ) *
                ( currentRotationMatrix.transpose( ) * currentThrust );
        dependentVariableResult[ outputIterator->first ] = newOutput;
    }

    // Write propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "thrustExampleFromFilePropagationHistory.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write dependent variable history to file.
    input_output::writeDataMapToTextFile( dependentVariableResult,
                                          "thrustExampleFromFileDependentVariableHistory.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );
}
