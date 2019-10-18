/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>
#include <Tudat/Basics/testMacros.h>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <SatellitePropagatorExamples/applicationOutput.h>
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/baseFunctionsSphericalShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/compositeFunctionSphericalShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/sphericalShaping.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"


int main( )
{
    using namespace tudat;
    using namespace tudat::input_output;
    using namespace tudat::simulation_setup;
    using namespace tudat::shape_based_methods;

    std::string outputSubFolder = "ShapeBasedTrajectoriesExample/";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////           DEFINE TRAJECTORY GLOBAL PARAMETERS      ////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int numberOfRevolutions = 1;
    double julianDateAtDeparture = 8174.5 * physical_constants::JULIAN_DAY;
    double  timeOfFlight = 580.0 * physical_constants::JULIAN_DAY;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Retrieve cartesian state at departure and arrival.
    Eigen::Vector6d cartesianStateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDateAtDeparture );
    Eigen::Vector6d cartesianStateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDateAtDeparture + timeOfFlight );



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////         DEFINE PERTURBED DYNAMICAL ENVIRONMENT          /////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    spice_interface::loadStandardSpiceKernels( );

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Borzi" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Borzi" )->setSuppressDependentOrientationCalculatorWarning( true );
    bodyMap.at( "Borzi" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Define body to propagate and central body.
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Borzi" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Sun" );

    // Set vehicle mass.
    bodyMap[ "Borzi" ]->setConstantBodyMass( 2000.0 );


    // Define specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return 3000.0;
    };

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Borzi" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Borzi", bodyMap ) );


    // Define acceleration map for the simplified problem.
    // (empty map as central gravity and low-thrust accelerations are already included in the shape-based methods)
    basic_astrodynamics::AccelerationMap perturbingAccelerationsMapSimplifiedProblem;


    // Define acceleration map for the fully perturbed problem.

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationSettingsPerturbedProblem;
    accelerationSettingsPerturbedProblem[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationSettingsPerturbedProblem[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationSettingsPerturbedProblem[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationSettingsPerturbedProblem[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );

    SelectedAccelerationMap accelerationMap;

    accelerationMap[ "Borzi" ] = accelerationSettingsPerturbedProblem;
    bodiesToPropagate.push_back( "Borzi" );
    centralBodies.push_back( "Sun" );

    basic_astrodynamics::AccelerationMap perturbingAccelerationsMapPertubedProblem = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////                 DEFINE PROPAGATION SETTINGS                   /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define integrator settings.
    double stepSize = timeOfFlight / static_cast< double >( 8000.0 );
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////             DEFINE HODOGRAPHIC SHAPING LEG          ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
    double scaleFactor = 1.0 / timeOfFlight;

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Create components for the axial velocity composite function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );


    // Initialize free coefficients vector for radial velocity function (empty here, only 3 base functions).
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for normal velocity function (empty here, only 3 base functions).
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for axial velocity function (empty here, only 3 base functions).
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    shape_based_methods::HodographicShaping hodographicShaping(
                cartesianStateAtDeparture, cartesianStateAtArrival, timeOfFlight, 1, bodyMap, "Borzi", "Sun",
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction,
                integratorSettings );


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////             DEFINE SPHERICAL SHAPING LEG            ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
            std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 30 );

    // Compute shaped trajectory.
    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                cartesianStateAtDeparture, cartesianStateAtArrival, timeOfFlight,
                numberOfRevolutions, bodyMap, "Borzi", "Sun", 0.000703,
                rootFinderSettings, 1.0e-6, 1.0e-1, integratorSettings );


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////       NUMERICALLY PROPAGATE THE SIMPLIFIED PROBLEM            /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map< double, Eigen::VectorXd > hodographicShapingPropagationUnperturbedCase;
    std::map< double, Eigen::Vector6d > hodographicShapingAnalyticalResults;
    std::map< double, Eigen::VectorXd > hodographicShapingDependentVariablesHistory;

    // Create propagator settings for hodographic shaping.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            hodographicShapingPropagatorSettings = hodographicShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapSimplifiedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    hodographicShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, hodographicShapingPropagatorSettings, hodographicShapingPropagationUnperturbedCase,
                hodographicShapingAnalyticalResults, hodographicShapingDependentVariablesHistory );


    std::map< double, Eigen::VectorXd > sphericalShapingPropagationUnperturbedCase;
    std::map< double, Eigen::Vector6d > sphericalShapingAnalyticalResults;
    std::map< double, Eigen::VectorXd > sphericalShapingDependentVariablesHistory;

    // Create propagator settings for spherical shaping.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            sphericalShapingPropagatorSettings = sphericalShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapSimplifiedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    sphericalShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, sphericalShapingPropagatorSettings, sphericalShapingPropagationUnperturbedCase,
                sphericalShapingAnalyticalResults, sphericalShapingDependentVariablesHistory );

    input_output::writeDataMapToTextFile( hodographicShapingAnalyticalResults,
                                          "hodographicShapingAnalyticalResults.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sphericalShapingAnalyticalResults,
                                          "sphericalShapingAnalyticalResults.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( hodographicShapingPropagationUnperturbedCase,
                                          "hodographicShapingPropagationUnperturbedCase.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sphericalShapingPropagationUnperturbedCase,
                                          "sphericalShapingPropagationUnperturbedCase.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////         RETRIEVE TRAJECTORY, MASS, THRUST AND THRUST ACCELERATION PROFILES           //////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Hodographic shaping
    std::vector< double > epochsVectorHodographicShaping;
    for ( std::map< double, Eigen::Vector6d >::iterator itr = hodographicShapingAnalyticalResults.begin( ) ;
          itr != hodographicShapingAnalyticalResults.end( ) ; itr++ )
    {
        epochsVectorHodographicShaping.push_back( itr->first );
    }

    std::map< double, Eigen::VectorXd > hodographicShapingMassProfile;
    std::map< double, Eigen::VectorXd > hodographicShapingThrustProfile;
    std::map< double, Eigen::VectorXd > hodographicShapingThrustAccelerationProfile;

    hodographicShaping.getMassProfile(
                epochsVectorHodographicShaping, hodographicShapingMassProfile, specificImpulseFunction, integratorSettings );
    hodographicShaping.getThrustProfile(
                epochsVectorHodographicShaping, hodographicShapingThrustProfile, specificImpulseFunction, integratorSettings );
    hodographicShaping.getThrustAccelerationProfile(
                epochsVectorHodographicShaping, hodographicShapingThrustAccelerationProfile, specificImpulseFunction, integratorSettings );


    // Spherical shaping
    std::vector< double > epochsVectorSphericalShaping;
    for ( std::map< double, Eigen::Vector6d >::iterator itr = sphericalShapingAnalyticalResults.begin( ) ;
          itr != sphericalShapingAnalyticalResults.end( ) ; itr++ )
    {
        epochsVectorSphericalShaping.push_back( itr->first );
    }

    std::map< double, Eigen::VectorXd > sphericalShapingMassProfile;
    std::map< double, Eigen::VectorXd > sphericalShapingThrustProfile;
    std::map< double, Eigen::VectorXd > sphericalShapingThrustAccelerationProfile;

    sphericalShaping.getMassProfile(
                epochsVectorSphericalShaping, sphericalShapingMassProfile, specificImpulseFunction, integratorSettings );
    sphericalShaping.getThrustProfile(
                epochsVectorSphericalShaping, sphericalShapingThrustProfile, specificImpulseFunction, integratorSettings );
    sphericalShaping.getThrustAccelerationProfile(
                epochsVectorSphericalShaping, sphericalShapingThrustAccelerationProfile, specificImpulseFunction, integratorSettings );



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////         PROPAGATE THE FULLY PERTURBED PROBLEM           /////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map< double, Eigen::VectorXd > hodographicShapingPropagationPerturbedCase;
    std::map< double, Eigen::Vector6d > hodographicShapingAnalyticalResultsPerturbedCase;
    std::map< double, Eigen::VectorXd > hodographicShapingDependentVariablesHistoryPerturbedCase;

    std::map< double, Eigen::VectorXd > sphericalShapingPropagationPerturbedCase;
    std::map< double, Eigen::Vector6d > sphericalShapingAnalyticalResultsPerturbedCase;
    std::map< double, Eigen::VectorXd > sphericalShapingDependentVariablesHistoryPerturbedCase;

    // Create propagator settings for hodographic shaping.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            hodographicShapingPropagatorSettingsPerturbedCase = hodographicShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapPertubedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    hodographicShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, hodographicShapingPropagatorSettingsPerturbedCase, hodographicShapingPropagationPerturbedCase,
                hodographicShapingAnalyticalResultsPerturbedCase, hodographicShapingDependentVariablesHistoryPerturbedCase );

    // Create propagator settings for spherical shaping.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            sphericalShapingPropagatorSettingsPerturbedCase = sphericalShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapPertubedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    sphericalShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, sphericalShapingPropagatorSettingsPerturbedCase, sphericalShapingPropagationPerturbedCase,
                sphericalShapingAnalyticalResultsPerturbedCase, sphericalShapingDependentVariablesHistoryPerturbedCase );



    input_output::writeDataMapToTextFile( hodographicShapingPropagationPerturbedCase,
                                          "hodographicShapingPropagationPerturbedCase.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sphericalShapingPropagationPerturbedCase,
                                          "sphericalShapingPropagationPerturbedCase.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( hodographicShapingMassProfile,
                                          "hodographicShapingMassProfile.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( hodographicShapingThrustProfile,
                                          "hodographicShapingThrustProfile.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( hodographicShapingThrustAccelerationProfile,
                                          "hodographicShapingThrustAccelerationProfile.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sphericalShapingMassProfile,
                                          "sphericalShapingMassProfile.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sphericalShapingThrustProfile,
                                          "sphericalShapingThrustProfile.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sphericalShapingThrustAccelerationProfile,
                                          "sphericalShapingThrustAccelerationProfile.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    std::cout << "deltaV hodographic shaping: " << hodographicShaping.computeDeltaV( ) << "\n\n";
    std::cout << "deltaV spherical shaping: " << sphericalShaping.computeDeltaV( ) << "\n\n";


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}