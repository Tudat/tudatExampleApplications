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
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/sphericalShaping.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"


int main( )
{
    using namespace tudat;
    using namespace tudat::input_output;
    using namespace tudat::simulation_setup;
    using namespace tudat::shape_based_methods;


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////           DEFINE TRAJECTORY GLOBAL PARAMETERS      ////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    double numberOfRevolutions = 1.0;
//    double julianDateAtDeparture = 2458849.5;
//    double timeOfFlight = 500.0 * physical_constants::JULIAN_DAY;

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

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );


    // Define body to propagate and central body.
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Vehicle" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Sun" );

//    // Acceleration from the central body.
//    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
//    bodyToPropagateAccelerations[ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
//                                                                basic_astrodynamics::central_gravity ) );

//    simulation_setup::SelectedAccelerationMap accelerationMap;
//    accelerationMap[ "Vehicle" ] = bodyToPropagateAccelerations;

//    // Create the acceleration map.
//    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
//                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Set vehicle mass.
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 2000.0 );


    // Define specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return 3000.0;
    };


    // Define acceleration map for the simplified problem.
    // (empty map as central gravity and low-thrust accelerations are already included in the shape-based methods)
    basic_astrodynamics::AccelerationMap perturbingAccelerationsMapSimplifiedProblem;



    // Define acceleration map for the perturbed problem.

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationSettingsPerturbedProblem;
    accelerationSettingsPerturbedProblem[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationSettingsPerturbedProblem[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationSettingsPerturbedProblem[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );

    SelectedAccelerationMap accelerationMap;

    accelerationMap[ "Vehicle" ] = accelerationSettingsPerturbedProblem;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Sun" );

    basic_astrodynamics::AccelerationMap perturbingAccelerationsMapPertubedProblem = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////                 DEFINE PROPAGATION SETTINGS                   /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define integrator settings.
    double stepSize = timeOfFlight / static_cast< double >( 20 );
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize / 400.0 );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_mass_rate_dependent_variables, "Vehicle" ) );

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
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, fourthNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, fifthNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                4.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                4.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, fourthAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, fifthAxialVelocityBaseFunctionSettings ) );

    // Initialize free coefficients vector for radial velocity function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsRadialVelocityFunction[ 0 ] = 500.0;
    freeCoefficientsRadialVelocityFunction[ 1 ] = 500.0;

    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsNormalVelocityFunction[ 0 ] = 500.0;
    freeCoefficientsNormalVelocityFunction[ 1 ] = - 200.0;

    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsAxialVelocityFunction[ 0 ] = 500.0;
    freeCoefficientsAxialVelocityFunction[ 1 ] = 2000.0;

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    shape_based_methods::HodographicShaping hodographicShaping(
                cartesianStateAtDeparture, cartesianStateAtArrival, timeOfFlight, 1, bodyMap, "Vehicle", "Sun",
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction,
                integratorSettings );

    std::cout << "hodographic shaping deltaV: " << hodographicShaping.computeDeltaV() << "\n\n";


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////             DEFINE SPHERICAL SHAPING LEG            ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 30 );

    // Compute shaped trajectory.
    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                cartesianStateAtDeparture, cartesianStateAtArrival, timeOfFlight,
                numberOfRevolutions, bodyMap, "Vehicle", "Sun", 0.000703,
                rootFinderSettings, 1.0e-6, 1.0e-1, integratorSettings );


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////       NUMERICALLY PROPAGATE THE SIMPLIFIED PROBLEM            /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map< double, Eigen::VectorXd > hodographicShapingFullPropagationResults;
    std::map< double, Eigen::Vector6d > hodographicShapingAnalyticalResults;
    std::map< double, Eigen::VectorXd > hodographicShapingDependentVariablesHistory;

    // Create propagator settings for hodographic shaping.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            hodographicShapingPropagatorSettings = hodographicShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapSimplifiedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    hodographicShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, hodographicShapingPropagatorSettings, hodographicShapingFullPropagationResults,
                hodographicShapingAnalyticalResults, hodographicShapingDependentVariablesHistory );


    std::map< double, Eigen::VectorXd > sphericalShapingFullPropagationResults;
    std::map< double, Eigen::Vector6d > sphericalShapingAnalyticalResults;
    std::map< double, Eigen::VectorXd > sphericalShapingDependentVariablesHistory;

    // Create propagator settings for spherical shaping.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            sphericalShapingPropagatorSettings = sphericalShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapSimplifiedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    sphericalShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, sphericalShapingPropagatorSettings, sphericalShapingFullPropagationResults,
                sphericalShapingAnalyticalResults, sphericalShapingDependentVariablesHistory );

    input_output::writeDataMapToTextFile( hodographicShapingAnalyticalResults,
                                          "hodographicShapingAnalyticalResults.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sphericalShapingAnalyticalResults,
                                          "sphericalShapingAnalyticalResults.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( hodographicShapingFullPropagationResults,
                                          "hodographicShapingFullPropagationResults.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sphericalShapingFullPropagationResults,
                                          "sphericalShapingFullPropagationResults.dat",
                                          tudat_applications::getOutputPath( ),
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

    std::map< double, Eigen::Vector6d > hodographicShapingTrajectory;
    std::map< double, Eigen::VectorXd > hodographicShapingMassProfile;
    std::map< double, Eigen::VectorXd > hodographicShapingThrustProfile;
    std::map< double, Eigen::VectorXd > hodographicShapingThrustAccelerationProfile;

//    hodographicShaping.getTrajectory( epochsVectorHodographicShaping, hodographicShapingTrajectory );
    hodographicShaping.getMassProfile(
                epochsVectorHodographicShaping, hodographicShapingMassProfile, specificImpulseFunction, integratorSettings );
//    hodographicShaping.getThrustProfile(
//                epochsVectorHodographicShaping, hodographicShapingThrustProfile, specificImpulseFunction, integratorSettings );
//    hodographicShaping.getThrustAccelerationProfile(
//                epochsVectorHodographicShaping, hodographicShapingThrustAccelerationProfile, specificImpulseFunction, integratorSettings );


    // Spherical shaping
    std::vector< double > epochsVectorSphericalShaping;
    for ( std::map< double, Eigen::Vector6d >::iterator itr = sphericalShapingAnalyticalResults.begin( ) ;
          itr != sphericalShapingAnalyticalResults.end( ) ; itr++ )
    {
        epochsVectorSphericalShaping.push_back( itr->first );
    }

    std::map< double, Eigen::Vector6d > sphericalShapingTrajectory;
    std::map< double, Eigen::VectorXd > sphericalShapingMassProfile;
    std::map< double, Eigen::VectorXd > sphericalShapingThrustProfile;
    std::map< double, Eigen::VectorXd > sphericalShapingThrustAccelerationProfile;

//    sphericalShaping.getTrajectory( epochsVectorSphericalShaping, sphericalShapingTrajectory );
    sphericalShaping.getMassProfile(
                epochsVectorSphericalShaping, sphericalShapingMassProfile, specificImpulseFunction, integratorSettings );
//    sphericalShaping.getThrustProfile(
//                epochsVectorSphericalShaping, sphericalShapingThrustProfile, specificImpulseFunction, integratorSettings );
//    sphericalShaping.getThrustAccelerationProfile(
//                epochsVectorSphericalShaping, sphericalShapingThrustAccelerationProfile, specificImpulseFunction, integratorSettings );



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////         PROPAGATE THE FULLY PERTURBED PROBLEM           /////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    hodographicShapingFullPropagationResults.clear( );
    hodographicShapingAnalyticalResults.clear( );
    hodographicShapingDependentVariablesHistory.clear( );

    // Create propagator settings for hodographic shaping.
//    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            hodographicShapingPropagatorSettings = hodographicShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapPertubedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    hodographicShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, hodographicShapingPropagatorSettings, hodographicShapingFullPropagationResults,
                hodographicShapingAnalyticalResults, hodographicShapingDependentVariablesHistory );


    sphericalShapingFullPropagationResults.clear( );
    sphericalShapingAnalyticalResults.clear( );
    sphericalShapingDependentVariablesHistory.clear( );

    // Create propagator settings for spherical shaping.
//    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            sphericalShapingPropagatorSettings = sphericalShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapPertubedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    sphericalShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, sphericalShapingPropagatorSettings, sphericalShapingFullPropagationResults,
                sphericalShapingAnalyticalResults, sphericalShapingDependentVariablesHistory );



    input_output::writeDataMapToTextFile( hodographicShapingFullPropagationResults,
                                          "hodographicShapingFullPropagationResultsPerturbedProblem.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( sphericalShapingFullPropagationResults,
                                          "sphericalShapingFullPropagationResultsPerturbedProblem.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    std::cout << "deltaV hodographic shaping: " << hodographicShaping.computeDeltaV( ) << "\n\n";
    std::cout << "deltaV spherical shaping: " << sphericalShaping.computeDeltaV( ) << "\n\n";






//    spice_interface::loadStandardSpiceKernels( );


//    // Global characteristics of the problem
//    double distanceSunJupiter = 778.0e9;

//    // Initialise the spacecraft state (B. Taylor, D. (1981). Horseshoe periodic orbits in the restricted problem of three bodies
//    // for a sun-Jupiter mass ratio. Astronomy and Astrophysics. 103. 288-294.)
//    Eigen::Vector6d initialState = Eigen::Vector6d::Zero();
//    initialState[0] = - 7.992e11;
//    initialState[4] =  -1.29e4;

//    // Create integrator settings.
//    double initialTime = 0.0;
//    const double fixedStepSize = 100000.0;
//    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
//            std::make_shared < numerical_integrators::IntegratorSettings < > >
//            ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize );

//    // Create body map.
//    std::vector < std::string > bodiesCR3BP;
//    bodiesCR3BP.push_back( "Sun" );
//    bodiesCR3BP.push_back( "Jupiter" );

//    // Define propagator settings variables.
//    std::vector< std::string > bodiesToPropagate;
//    std::vector< std::string > centralBodies;
//    bodiesToPropagate.push_back( "Spacecraft" );
//    centralBodies.push_back( "SSB" );

//    // Define final time for the propagation.
//    double gravitationalParameterSun = createGravityFieldModel(
//                getDefaultGravityFieldSettings(
//                    "Sun", TUDAT_NAN, TUDAT_NAN ), "Sun" )->getGravitationalParameter( );
//    double gravitationalParameterJupiter = createGravityFieldModel(
//                getDefaultGravityFieldSettings(
//                    "Jupiter", TUDAT_NAN, TUDAT_NAN ), "Jupiter" )->getGravitationalParameter( );
//    double finalTime = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(
//                29.2386 * ( 2.0 * mathematical_constants::PI ), gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter);

//    NamedBodyMap idealBodyMap = propagators::setupBodyMapCR3BP(
//                distanceSunJupiter, "Sun", "Jupiter", "Spacecraft" );

//    std::map< double, Eigen::Vector6d> fullPropagation;
//    std::map< double, Eigen::Vector6d> cr3bpPropagation;

//    /// Ideal case: full dynamics problem with the CR3BP assumtions
//    {


//        // Create acceleration map.
//        basic_astrodynamics::AccelerationMap accelerationModelMap = propagators::setupAccelerationMapCR3BP(
//                    "Sun", "Jupiter", bodiesToPropagate.at( 0 ), centralBodies.at( 0 ), idealBodyMap );

//        // Calculate the difference between CR3BP and full problem.
//        propagators::propagateCR3BPAndFullDynamicsProblem(
//                    initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
//                    bodiesToPropagate, centralBodies, idealBodyMap, bodiesCR3BP, fullPropagation,
//                    cr3bpPropagation );

//        Eigen::Vector6d stateDifference =
//                fullPropagation.rbegin( )->second - cr3bpPropagation.rbegin( )->second;

//        std::cout << "state difference at final time: " << stateDifference << std::endl;
//    }


//    std::map< double, Eigen::Vector6d> fullPropagationPerturbedCase;
//    std::map< double, Eigen::Vector6d> cr3bpPropagationPerturbedCase;

//    /// Perturbed case
//    {

//        std::string frameOrigin = "SSB";
//        std::string frameOrientation = "ECLIPJ2000";

//        NamedBodyMap perturbedBodyMap;


//        std::vector< std::string > additionalBodies = { "Earth", "Mars", "Venus", "Saturn" };
//        for( unsigned int i = 0; i < additionalBodies.size( ); i++ )
//        {

//            perturbedBodyMap[ additionalBodies.at( i ) ] = std::make_shared< Body >( );
//            perturbedBodyMap[ additionalBodies.at( i ) ]->setEphemeris(
//                        std::make_shared< ephemerides::ApproximatePlanetPositions>(
//                                additionalBodies.at( i ) ) );
//            perturbedBodyMap[ additionalBodies.at( i ) ]->setGravityFieldModel(
//                        createGravityFieldModel(
//                            std::make_shared< CentralGravityFieldSettings >(
//                                spice_interface::getBodyGravitationalParameter(
//                                    additionalBodies.at( i ) ) ), additionalBodies.at( i ) ) );
//        }

//        perturbedBodyMap[ "Sun" ] = idealBodyMap[ "Sun" ];
//        perturbedBodyMap[ "Jupiter" ] = idealBodyMap[ "Jupiter" ];


//        // Create the body to be propagated.
//        perturbedBodyMap[ "Spacecraft" ] = std::make_shared< Body >( );
//        perturbedBodyMap[ "Spacecraft" ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
//                                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                                            < double, Eigen::Vector6d > >( ), "SSB", frameOrientation ) );

//        setGlobalFrameBodyEphemerides( perturbedBodyMap, frameOrigin, frameOrientation );


//        // Set of accelerations experienced by the spacecraft.
//        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > bodyToPropagateAccelerations;
//        bodyToPropagateAccelerations["Sun"].push_back(std::make_shared< AccelerationSettings >(
//                                                          basic_astrodynamics::central_gravity ) );
//        bodyToPropagateAccelerations["Jupiter"].push_back(std::make_shared< AccelerationSettings >(
//                                                              basic_astrodynamics::central_gravity ) );
//        bodyToPropagateAccelerations["Earth"].push_back(std::make_shared< AccelerationSettings >(
//                                                            basic_astrodynamics::central_gravity ) );
//        bodyToPropagateAccelerations["Mars"].push_back(std::make_shared< AccelerationSettings >(
//                                                           basic_astrodynamics::central_gravity ) );
//        bodyToPropagateAccelerations["Venus"].push_back(std::make_shared< AccelerationSettings >(
//                                                            basic_astrodynamics::central_gravity ) );
//        bodyToPropagateAccelerations["Saturn"].push_back(std::make_shared< AccelerationSettings >(
//                                                             basic_astrodynamics::central_gravity ) );

//        SelectedAccelerationMap accelerationMap;
//        accelerationMap[ "Spacecraft" ] = bodyToPropagateAccelerations;


//        // Create the acceleration map.
//        basic_astrodynamics::AccelerationMap accelerationModelMapPerturbedCase = createAccelerationModelsMap(
//                    perturbedBodyMap, accelerationMap, bodiesToPropagate, centralBodies );


//        // Calculate the difference between CR3BP and full problem.
//        propagators::propagateCR3BPAndFullDynamicsProblem(
//                    initialTime, finalTime, initialState, integratorSettings,
//                    accelerationModelMapPerturbedCase,
//                    bodiesToPropagate, centralBodies, perturbedBodyMap, bodiesCR3BP,
//                    fullPropagationPerturbedCase,
//                    cr3bpPropagationPerturbedCase );

//        Eigen::Vector6d stateDifferencePerturbedCase =
//                fullPropagationPerturbedCase.rbegin( )->second - cr3bpPropagationPerturbedCase.rbegin( )->second;

//        std::cout << "state difference at final time for the perturbed case: " << stateDifferencePerturbedCase << std::endl;

//    }

//    /// Ouputs

//    // Outputs for the ideal case
//    {
//        std::map< double, Eigen::Vector6d > fullPropagationNormalisedCoRotatingFrame;
//        for( std::map< double, Eigen::Vector6d >::iterator itr = fullPropagation.begin( );
//             itr != fullPropagation.end( ); itr++ ){
//            fullPropagationNormalisedCoRotatingFrame[ itr->first ] = tudat::circular_restricted_three_body_problem::convertCartesianToCorotatingNormalizedCoordinates(
//                        gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter, itr->second, itr->first);
//        }

//        std::map< double, Eigen::Vector6d > cr3bpNormalisedCoRotatingFrame;
//        for( std::map< double, Eigen::Vector6d >::iterator itr = cr3bpPropagation.begin( );
//             itr != cr3bpPropagation.end( ); itr++ ){
//            cr3bpNormalisedCoRotatingFrame[ itr->first ] = tudat::circular_restricted_three_body_problem::convertCartesianToCorotatingNormalizedCoordinates(
//                        gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter, itr->second, itr->first);
//        }


//        input_output::writeDataMapToTextFile( fullPropagation,
//                                              "fullProblemPropagation.dat",
//                                              tudat_applications::getOutputPath( ),
//                                              "",
//                                              std::numeric_limits< double >::digits10,
//                                              std::numeric_limits< double >::digits10,
//                                              "," );

//        input_output::writeDataMapToTextFile( fullPropagationNormalisedCoRotatingFrame,
//                                              "fullProblemPropagationNormalisedCoRotatingFrame.dat",
//                                              tudat_applications::getOutputPath( ),
//                                              "",
//                                              std::numeric_limits< double >::digits10,
//                                              std::numeric_limits< double >::digits10,
//                                              "," );

//        input_output::writeDataMapToTextFile( cr3bpPropagation,
//                                              "CR3BPsolution.dat",
//                                              tudat_applications::getOutputPath( ),
//                                              "",
//                                              std::numeric_limits< double >::digits10,
//                                              std::numeric_limits< double >::digits10,
//                                              "," );

//        input_output::writeDataMapToTextFile( cr3bpNormalisedCoRotatingFrame,
//                                              "CR3BPnormalisedCoRotatingFrame.dat",
//                                              tudat_applications::getOutputPath( ),
//                                              "",
//                                              std::numeric_limits< double >::digits10,
//                                              std::numeric_limits< double >::digits10,
//                                              "," );
//    }

//    // Outputs for the perturbed case
//    {
//        std::map< double, Eigen::Vector6d > fullPropagationNormalisedCoRotatingFramePerturbedCase;
//        for( std::map< double, Eigen::Vector6d >::iterator itr = fullPropagationPerturbedCase.begin( );
//             itr != fullPropagationPerturbedCase.end( ); itr++ ){
//            fullPropagationNormalisedCoRotatingFramePerturbedCase[ itr->first ] = tudat::circular_restricted_three_body_problem::
//                    convertCartesianToCorotatingNormalizedCoordinates(gravitationalParameterSun, gravitationalParameterJupiter,
//                                                                      distanceSunJupiter, itr->second, itr->first);
//        }

//        std::map< double, Eigen::Vector6d > cr3bpNormalisedCoRotatingFramePerturbedCase;
//        for( std::map< double, Eigen::Vector6d >::iterator itr = cr3bpPropagationPerturbedCase.begin( );
//             itr != cr3bpPropagationPerturbedCase.end( ); itr++ ){
//            cr3bpNormalisedCoRotatingFramePerturbedCase[ itr->first ] = tudat::circular_restricted_three_body_problem::
//                    convertCartesianToCorotatingNormalizedCoordinates(
//                        gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter, itr->second, itr->first);
//        }


//        input_output::writeDataMapToTextFile( fullPropagationPerturbedCase,
//                                              "fullProblemPropagationPerturbedCase.dat",
//                                              tudat_applications::getOutputPath( ),
//                                              "",
//                                              std::numeric_limits< double >::digits10,
//                                              std::numeric_limits< double >::digits10,
//                                              "," );

//        input_output::writeDataMapToTextFile( fullPropagationNormalisedCoRotatingFramePerturbedCase,
//                                              "fullProblemPropagationNormalisedCoRotatingFramePerturbedCase.dat",
//                                              tudat_applications::getOutputPath( ),
//                                              "",
//                                              std::numeric_limits< double >::digits10,
//                                              std::numeric_limits< double >::digits10,
//                                              "," );

//        input_output::writeDataMapToTextFile( cr3bpPropagationPerturbedCase,
//                                              "CR3BPsolutionPerturbedCase.dat",
//                                              tudat_applications::getOutputPath( ),
//                                              "",
//                                              std::numeric_limits< double >::digits10,
//                                              std::numeric_limits< double >::digits10,
//                                              "," );

//        input_output::writeDataMapToTextFile( cr3bpNormalisedCoRotatingFramePerturbedCase,
//                                              "CR3BPnormalisedCoRotatingFramePerturbedCase.dat",
//                                              tudat_applications::getOutputPath( ),
//                                              "",
//                                              std::numeric_limits< double >::digits10,
//                                              std::numeric_limits< double >::digits10,
//                                              "," );
//    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
