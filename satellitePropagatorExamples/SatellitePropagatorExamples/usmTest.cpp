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
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "matlabtest.cpp" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

//! Execute propagation of orbit of Satellite around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::ephemerides;
    using namespace tudat::input_output;
    using namespace tudat::aerodynamics;
    using namespace tudat::root_finders;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;//7.0 * tudat::physical_constants::JULIAN_YEAR + 30.0 * 6.0 * tudat::physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = 1.0 * tudat::physical_constants::JULIAN_DAY + simulationStartEpoch;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Earth" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    bodySettings[ "Mars" ]->gravityFieldSettings = std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    std::string atmosphereFile =  input_output::getTudatRootPath( ) + "/External/AtmosphereTables/MCDMeanAtmosphere.dat";
    bodySettings[ "Mars" ]->atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >( atmosphereFile );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Satellite" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Satellite" ]->setConstantBodyMass( 300.0 );

    // Aerodynamic parameters
    const double referenceAreaAerodynamic = 37.5;
    Eigen::Vector3d aerodynamicCoefficient = Eigen::Vector3d::Zero( );
    aerodynamicCoefficient( 0 ) = 2.2; // drag coefficient

    // Create aerodynamic coefficient interface settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceAreaAerodynamic, aerodynamicCoefficient, 1, 1 );

    bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 20;
    double radiationPressureCoefficient = 1.25;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    std::shared_ptr< RadiationPressureInterfaceSettings > SatelliteRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Satellite" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    SatelliteRadiationPressureSettings, "Satellite", bodyMap ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    accelerationsOfSatellite[ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 21, 21 ) );

    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        if( i != 1 )
        {
            accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( std::make_shared< AccelerationSettings >(
                                                                              basic_astrodynamics::central_gravity ) );
        }
    }
    accelerationsOfSatellite[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );

    accelerationsOfSatellite[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::aerodynamic ) );

    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    bodiesToPropagate.push_back( "Satellite" );
    centralBodies.push_back( "Mars" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Satellite.
    Eigen::Vector6d SatelliteInitialStateInKeplerianElements;
    SatelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 23477500;//25945000;
    SatelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.848152;//0.865099;
    SatelliteInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 93.0 );
    SatelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 158.7 );
    SatelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    SatelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 180.0 );

    double earthGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d SatelliteInitialState = convertKeplerianToCartesianElements(
                SatelliteInitialStateInKeplerianElements, earthGravitationalParameter );

    // Propagation termination conditions.
    std::shared_ptr< PropagationTerminationSettings > terminationSettingsTime = std::make_shared<
            PropagationTimeTerminationSettings >( simulationEndEpoch );
    std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettingAltitude =
            std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Satellite", "Mars" );
    std::shared_ptr< PropagationTerminationSettings > terminationSettingsAltitude = std::make_shared<
            PropagationDependentVariableTerminationSettings >( dependentVariableSettingAltitude, 10.0E+03, true );
    std::vector< std::shared_ptr< PropagationTerminationSettings > > propagationTerminationSettings;
    propagationTerminationSettings.push_back( terminationSettingsTime );
    propagationTerminationSettings.push_back( terminationSettingsAltitude );

    std::shared_ptr< PropagationTerminationSettings > terminationSettings =
            std::make_shared< PropagationHybridTerminationSettings >( propagationTerminationSettings, true );

    // Dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable, "Satellite", "Mars" ) );
    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable, "Satellite", "Mars" ) );
    dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::aerodynamic, "Satellite", "Mars", false ) );

    // Propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, SatelliteInitialState, terminationSettings,
              cowell, std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

    // Integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, 1.0 );
//    std::shared_ptr< IntegratorSettings< > > integratorSettings =
//            std::make_shared< RungeKuttaVariableStepSizeSettings< > >(
//                rungeKuttaVariableStepSize, simulationStartEpoch, 1.0,
//                RungeKuttaCoefficients::rungeKuttaFehlberg56, 0.005, 10.0, 1e-15, 1e-15 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > keplerianIntegrationResult;
    std::map< double, Eigen::VectorXd > dependentVariableSolution = dynamicsSimulator.getDependentVariableHistory( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Compute map of Kepler elements
    Eigen::Vector6d currentCartesianState;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
         stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
    {
        // Retrieve current Cartesian state (convert to Earth-centered frame if needed)
        currentCartesianState = stateIterator->second;

        if( centralBodies.at( 0 ) == "SSB" )
        {
            currentCartesianState -=
                    bodyMap[ "Mars" ]->getStateInBaseFrameFromEphemeris( stateIterator->first );
        }

        keplerianIntegrationResult[ stateIterator->first ] =
                convertCartesianToKeplerianElements( currentCartesianState, earthGravitationalParameter );
    }

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( cartesianIntegrationResult,
                                          "trajectory.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( keplerianIntegrationResult,
                                          "orbit.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( dependentVariableSolution,
                                          "dependent.dat",getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
