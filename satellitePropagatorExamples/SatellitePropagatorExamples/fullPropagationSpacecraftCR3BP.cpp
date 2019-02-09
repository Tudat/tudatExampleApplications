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

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Basics/testMacros.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationCR3BPFullProblem.h"

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <SatellitePropagatorExamples/applicationOutput.h>
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"


int main( )
{
    using namespace tudat;
    using namespace tudat::input_output;



    // Global characteristics of the problem
    double distanceSunJupiter = 778.0e9;

    // Initialise the spacecraft state (B. Taylor, D. (1981). Horseshoe periodic orbits in the restricted problem of three bodies
    //                                 for a sun-Jupiter mass ratio. Astronomy and Astrophysics. 103. 288-294.)
    Eigen::Vector6d initialState = Eigen::Vector6d::Zero();
    initialState[0] = - 7.992e11;
    initialState[4] =  -1.29e4;


    // Create integrator settings.
    double initialTime = 0.0;
    const double fixedStepSize = 100000.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
            ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize );



    /// Ideal case: full dynamics problem with the CR3BP assumtions

    // Create body map.
    std::vector < std::string > bodiesCR3BP;
    bodiesCR3BP.push_back( "Sun" );
    bodiesCR3BP.push_back( "Jupiter" );

    simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapCR3BP(
                distanceSunJupiter, "Sun", "Jupiter", "Spacecraft" );


    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Spacecraft" );
    centralBodies.push_back( "SSB" );

    // Create acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = propagators::setupAccelerationMapCR3BP(
                "Sun", "Jupiter", bodiesToPropagate.at( 0 ), centralBodies.at( 0 ), bodyMap );

    // Define final time for the propagation.
    double gravitationalParameterSun = simulation_setup::createGravityFieldModel( simulation_setup::getDefaultGravityFieldSettings(
                    "Sun", TUDAT_NAN, TUDAT_NAN ), "Sun" )->getGravitationalParameter( );
    double gravitationalParameterJupiter = simulation_setup::createGravityFieldModel( simulation_setup::getDefaultGravityFieldSettings(
                    "Jupiter", TUDAT_NAN, TUDAT_NAN ), "Jupiter" )->getGravitationalParameter( );
    double finalTime = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(
                29.2386 * (2*mathematical_constants::PI), gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter);


    // Calculate the difference between CR3BP and full problem.
    std::map< double, Eigen::Vector6d> fullPropagation;
    std::map< double, Eigen::Vector6d> cr3bpPropagation;

    propagators::propagateCR3BPAndFullDynamicsProblem(initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
                                                      bodiesToPropagate, centralBodies, bodyMap, bodiesCR3BP, fullPropagation,
                                                      cr3bpPropagation);

    Eigen::Vector6d stateDifference = propagators::getFinalStateDifferenceFullPropagationWrtCR3BP(
                initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
                bodiesToPropagate, centralBodies,bodyMap, bodiesCR3BP );


    std::cout << "state difference at final time: " << stateDifference << std::endl;



    /// Perturbed case

    spice_interface::loadStandardSpiceKernels( );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "J2000";


    // Create body map.
    simulation_setup::NamedBodyMap bodyMapPerturbedCase;

    bodyMapPerturbedCase[ "Jupiter" ] = std::make_shared< simulation_setup::Body >( );
    bodyMapPerturbedCase[ "Jupiter" ]->setEphemeris( std::make_shared< ephemerides::ApproximatePlanetPositions>(
                        ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter) );
    bodyMapPerturbedCase[ "Jupiter" ]->setGravityFieldModel( simulation_setup::createGravityFieldModel(
                  simulation_setup::getDefaultGravityFieldSettings("Jupiter", TUDAT_NAN, TUDAT_NAN), "Jupiter"));


    // Ensure that the Sun and Jupiter are initially aligned.
    double distanceJupiterBarycenterPerturbedCase = bodyMapPerturbedCase[ "Jupiter" ]->getEphemeris()->getCartesianState(initialTime)
            .segment(0,3).norm();
    double distanceBarycenterSunPerturbedCase = distanceSunJupiter - distanceJupiterBarycenterPerturbedCase;

    Eigen::Vector3d initialPositionSun;
    initialPositionSun[0] = -( distanceBarycenterSunPerturbedCase / distanceJupiterBarycenterPerturbedCase ) * bodyMapPerturbedCase[ "Jupiter" ]
            ->getEphemeris()->getCartesianState(initialTime)[0];
    initialPositionSun[1] = -( distanceBarycenterSunPerturbedCase / distanceJupiterBarycenterPerturbedCase ) * bodyMapPerturbedCase[ "Jupiter" ]
            ->getEphemeris()->getCartesianState(initialTime)[1];
    initialPositionSun[2] = bodyMapPerturbedCase[ "Jupiter" ]->getEphemeris()->getCartesianState(initialTime)[2];

    // Define ephemeris for the Sun.
    Eigen::Vector6d initialStateSun = Eigen::Vector6d::Zero();
    initialStateSun.segment(0,3) = initialPositionSun;

    bodyMapPerturbedCase[ "Sun" ] = std::make_shared< simulation_setup::Body >( );
    bodyMapPerturbedCase[ "Sun" ]->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris> (
                                                     initialStateSun, frameOrigin, frameOrientation ) );
    bodyMapPerturbedCase[ "Sun" ]->setGravityFieldModel( simulation_setup::createGravityFieldModel(
                  simulation_setup::getDefaultGravityFieldSettings("Sun", TUDAT_NAN, TUDAT_NAN), "Sun"));



    // Create the body to be propagated.
    bodyMapPerturbedCase[ "Spacecraft" ] = std::make_shared< simulation_setup::Body >( );
    bodyMapPerturbedCase[ "Spacecraft" ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                      std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                      < double, Eigen::Vector6d > >( ), "SSB", frameOrientation ) );

    setGlobalFrameBodyEphemerides( bodyMapPerturbedCase, "SSB", frameOrientation );


    // Set of accelerations experienced by the spacecraft.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations["Sun"].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations["Jupiter"].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Spacecraft" ] = bodyToPropagateAccelerations;


    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMapPerturbedCase = createAccelerationModelsMap(
                bodyMapPerturbedCase, accelerationMap, bodiesToPropagate, centralBodies );


    // Calculate the difference between CR3BP and full problem.
    std::map< double, Eigen::Vector6d> fullPropagationPerturbedCase;
    std::map< double, Eigen::Vector6d> cr3bpPropagationPerturbedCase;

    propagators::propagateCR3BPAndFullDynamicsProblem(initialTime, finalTime, initialState, integratorSettings,
                                                      accelerationModelMapPerturbedCase,
                                                      bodiesToPropagate, centralBodies, bodyMapPerturbedCase, bodiesCR3BP,
                                                      fullPropagationPerturbedCase,
                                                      cr3bpPropagationPerturbedCase);

    Eigen::Vector6d stateDifferencePerturbedCase = propagators::getFinalStateDifferenceFullPropagationWrtCR3BP(
                initialTime, finalTime, initialState, integratorSettings, accelerationModelMapPerturbedCase,
                bodiesToPropagate, centralBodies, bodyMapPerturbedCase, bodiesCR3BP );

    std::cout << "state difference at final time for the perturbed case: " << stateDifferencePerturbedCase << std::endl;



    /// Ouputs

    // Outputs for the ideal case

    std::map< double, Eigen::Vector6d > fullPropagationNormalisedCoRotatingFrame;
    for( std::map< double, Eigen::Vector6d >::iterator itr = fullPropagation.begin( );
         itr != fullPropagation.end( ); itr++ ){
        fullPropagationNormalisedCoRotatingFrame[ itr->first ] = tudat::circular_restricted_three_body_problem::convertCartesianToCorotatingNormalizedCoordinates(
                    gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter, itr->second, itr->first);
    }

    std::map< double, Eigen::Vector6d > cr3bpNormalisedCoRotatingFrame;
    for( std::map< double, Eigen::Vector6d >::iterator itr = cr3bpPropagation.begin( );
         itr != cr3bpPropagation.end( ); itr++ ){
        cr3bpNormalisedCoRotatingFrame[ itr->first ] = tudat::circular_restricted_three_body_problem::convertCartesianToCorotatingNormalizedCoordinates(
                    gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter, itr->second, itr->first);
    }


    input_output::writeDataMapToTextFile( fullPropagation,
                                          "fullProblemPropagation.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( fullPropagationNormalisedCoRotatingFrame,
                                          "fullProblemPropagationNormalisedCoRotatingFrame.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( cr3bpPropagation,
                                          "CR3BPsolution.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( cr3bpNormalisedCoRotatingFrame,
                                          "CR3BPnormalisedCoRotatingFrame.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Outputs for the perturbed case

    std::map< double, Eigen::Vector6d > fullPropagationNormalisedCoRotatingFramePerturbedCase;
    for( std::map< double, Eigen::Vector6d >::iterator itr = fullPropagationPerturbedCase.begin( );
         itr != fullPropagationPerturbedCase.end( ); itr++ ){
        fullPropagationNormalisedCoRotatingFramePerturbedCase[ itr->first ] = tudat::circular_restricted_three_body_problem::
                convertCartesianToCorotatingNormalizedCoordinates(gravitationalParameterSun, gravitationalParameterJupiter,
                                                                  distanceSunJupiter, itr->second, itr->first);
    }

    std::map< double, Eigen::Vector6d > cr3bpNormalisedCoRotatingFramePerturbedCase;
    for( std::map< double, Eigen::Vector6d >::iterator itr = cr3bpPropagationPerturbedCase.begin( );
         itr != cr3bpPropagationPerturbedCase.end( ); itr++ ){
        cr3bpNormalisedCoRotatingFramePerturbedCase[ itr->first ] = tudat::circular_restricted_three_body_problem::
                convertCartesianToCorotatingNormalizedCoordinates(
                    gravitationalParameterSun, gravitationalParameterJupiter, distanceSunJupiter, itr->second, itr->first);
    }


    input_output::writeDataMapToTextFile( fullPropagationPerturbedCase,
                                          "fullProblemPropagationPerturbedCase.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( fullPropagationNormalisedCoRotatingFramePerturbedCase,
                                          "fullProblemPropagationNormalisedCoRotatingFramePerturbedCase.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( cr3bpPropagationPerturbedCase,
                                          "CR3BPsolutionPerturbedCase.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( cr3bpNormalisedCoRotatingFramePerturbedCase,
                                          "CR3BPnormalisedCoRotatingFramePerturbedCase.dat",
                                          tudat_applications::getOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
