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

#include <Tudat/InputOutput/basicInputOutput.h>
#include <SatellitePropagatorExamples/applicationOutput.h>

#include "Tudat/SimulationSetup/PropagationSetup/propagationPatchedConicFullProblem.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationLambertTargeterFullProblem.h"


#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/planetTrajectory.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

int main( )
{
    using namespace tudat;
    using namespace tudat::input_output;
    using namespace tudat::input_output::parsed_data_vector_utilities;
    using namespace tudat::transfer_trajectories;


    ///  Characteristics of the interplanetary trajectory

    // Define central body of the trajectory and body to be propagated.
    std::vector< std::string > centralBody; centralBody.push_back( "Sun" );
    std::string bodyToPropagate = "spacecraft";

    // Specify the number and types of legs and type of legs.
    int numberOfLegs = 5;
    std::vector< TransferLegType > legTypeVector( numberOfLegs );
    legTypeVector[ 0 ] = mga1DsmVelocity_Departure;
    legTypeVector[ 1 ] = mga1DsmVelocity_Swingby;
    legTypeVector[ 2 ] = mga1DsmVelocity_Swingby;
    legTypeVector[ 3 ] = mga1DsmVelocity_Swingby;
    legTypeVector[ 4 ] = capture;


    // Name of the bodies involved in the trajectory
    std::vector< std::string > transferBodyTrajectory;
    transferBodyTrajectory.push_back("Earth");
    transferBodyTrajectory.push_back("Earth");
    transferBodyTrajectory.push_back("Venus");
    transferBodyTrajectory.push_back("Venus");
    transferBodyTrajectory.push_back("Mercury");


    // Create variable vector.
    std::vector< double > variableVector;

    // Add the time of flight and start epoch.
    variableVector.push_back( 1171.64503236 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 399.999999715 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 178.372255301 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 299.223139512 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 180.510754824 * physical_constants::JULIAN_DAY);
    variableVector.push_back( 1.0); // The capture time is irrelevant for the final leg.

    // Add the additional variables.
    // 1st leg.
    variableVector.push_back( 0.234594654679 );
    variableVector.push_back( 1408.99421278 );
    variableVector.push_back( 0.37992647165 * 2 * 3.14159265358979 );
    variableVector.push_back( std::acos(  2 * 0.498004040298 - 1. ) - 3.14159265358979 / 2 );
    // 2nd leg.
    variableVector.push_back( 0.0964769387134 );
    variableVector.push_back( 1.35077257078 );
    variableVector.push_back( 1.80629232251 * 6.378e6 );
    variableVector.push_back( 0.0 );
    // 3rd leg.
    variableVector.push_back( 0.829948744508);
    variableVector.push_back( 1.09554368115 );
    variableVector.push_back( 3.04129845698 * 6.052e6 );
    variableVector.push_back( 0.0 );
    // 4th leg.
    variableVector.push_back( 0.317174785637 );
    variableVector.push_back( 1.34317576594 );
    variableVector.push_back( 1.10000000891 * 6.052e6 );
    variableVector.push_back( 0.0 );


    // Create minimum pericenter radii vector
    std::vector< double > minimumPericenterRadii;
    minimumPericenterRadii.push_back( TUDAT_NAN ); minimumPericenterRadii.push_back( TUDAT_NAN ); minimumPericenterRadii.push_back( TUDAT_NAN );
    minimumPericenterRadii.push_back( TUDAT_NAN ); minimumPericenterRadii.push_back( TUDAT_NAN );

    // Create departure and capture variables.
    std::vector< double > semiMajorAxes;
    semiMajorAxes.push_back( std::numeric_limits< double >::infinity( ) ); semiMajorAxes.push_back( std::numeric_limits< double >::infinity( ) );
    std::vector< double > eccentricities;
    eccentricities.push_back( 0.0 ); eccentricities.push_back( 0.0 );

    // Define integrator settings.
    double initialTime = 0.0;
    double fixedStepSize = 100.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > > ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize);



    /// Ideal case

    // Create body map.
    std::vector< double > gravitationalParametersTransferBodies;
    gravitationalParametersTransferBodies.push_back( 3.9860119e14 );
    gravitationalParametersTransferBodies.push_back( 3.9860119e14 );
    gravitationalParametersTransferBodies.push_back( 3.24860e14 );
    gravitationalParametersTransferBodies.push_back( 3.24860e14 );
    gravitationalParametersTransferBodies.push_back( 2.2321e13 );

    std::vector< ephemerides::EphemerisPointer > ephemerisVectorTransferBodies;
    ephemerisVectorTransferBodies.push_back( std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter) ) ;
    ephemerisVectorTransferBodies.push_back( std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter) ) ;
    ephemerisVectorTransferBodies.push_back( std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus) );
    ephemerisVectorTransferBodies.push_back( std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus) );
    ephemerisVectorTransferBodies.push_back( std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury) );

    simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapFromUserDefinedEphemeridesForPatchedConicsTrajectory(
            centralBody[0], bodyToPropagate, transferBodyTrajectory,
            ephemerisVectorTransferBodies, gravitationalParametersTransferBodies);


    // Create acceleration map.
    std::vector< basic_astrodynamics::AccelerationMap > accelerationMap = propagators::setupAccelerationMapPatchedConicsTrajectory(
                transferBodyTrajectory.size(), centralBody[0], bodyToPropagate, bodyMap);

    // Calculate the patched conics solution and the propagation results of the associated full dynamics problem for each leg.
    std::map< int, std::map< double, Eigen::Vector6d > > patchedConicsTrajectory;
    std::map< int, std::map< double, Eigen::Vector6d > > fullProblemTrajectory;

    propagators::fullPropagationPatchedConicsTrajectory( bodyMap, accelerationMap, transferBodyTrajectory, centralBody[0], bodyToPropagate,
            legTypeVector, variableVector, minimumPericenterRadii, semiMajorAxes, eccentricities, integratorSettings,
            patchedConicsTrajectory, fullProblemTrajectory, false);


    // Compute difference between patched conics trajectory and full problem at departure and at arrival for each leg.
    std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > differenceStateArrivalAndDeparturePerLeg =
            propagators::getDifferenceFullProblemWrtPatchedConicsTrajectory( bodyMap, accelerationMap, transferBodyTrajectory,
                centralBody[0], bodyToPropagate, legTypeVector, variableVector, minimumPericenterRadii, semiMajorAxes, eccentricities,
                integratorSettings, false );

    for( std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > >::iterator itr = differenceStateArrivalAndDeparturePerLeg.begin( );
         itr != differenceStateArrivalAndDeparturePerLeg.end( ); itr++ ){
        std::cout << "Leg " << itr->first << "\n\n";
        std::cout << "Departure: " << itr->second.first << "\n\n";
        std::cout << "Arrival: " << itr->second.second << "\n\n";
    }



    /// Perturbed case

    // Define accelerations acting on the spacecraft.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Sun" ].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ "Earth" ].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ "Venus" ].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ "Earth" ].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ "Sun" ].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMapPerturbedCase;
    accelerationMapPerturbedCase[ "spacecraft" ] = bodyToPropagateAccelerations;


    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMapPerturbedCase = simulation_setup::createAccelerationModelsMap(
                bodyMap, accelerationMapPerturbedCase, {bodyToPropagate}, {centralBody} );
    std::vector<  basic_astrodynamics::AccelerationMap > accelerationMapVectorPerturbedCase;
    for (int i = 0 ; i < numberOfLegs ; i++){
        accelerationMapVectorPerturbedCase.push_back( accelerationModelMapPerturbedCase );
    }


    // Calculate the patched conics trajectory and propagate the full dynamics problem jointly.
    std::map< int, std::map< double, Eigen::Vector6d > > patchedConicsTrajectoryPerturbedCase;
    std::map< int, std::map< double, Eigen::Vector6d > > fullProblemTrajectoryPerturbedCase;

    propagators::fullPropagationPatchedConicsTrajectory( bodyMap, accelerationMapVectorPerturbedCase,
            transferBodyTrajectory, centralBody[0], bodyToPropagate, legTypeVector, variableVector, minimumPericenterRadii,
            semiMajorAxes, eccentricities, integratorSettings, patchedConicsTrajectoryPerturbedCase,
            fullProblemTrajectoryPerturbedCase, false);

    // Compute difference between patched conics trajectory and full problem at departure and at arrival for each leg.
    std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > differenceStatePerturbedCase =
            propagators::getDifferenceFullProblemWrtPatchedConicsTrajectory( bodyMap, accelerationModelMapPerturbedCase, transferBodyTrajectory,
                centralBody[0], bodyToPropagate, legTypeVector, variableVector, minimumPericenterRadii, semiMajorAxes, eccentricities,
                integratorSettings, false );

    std::cout << "State difference perturbed case: " << "\n\n";
    for( std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > >::iterator itr = differenceStatePerturbedCase.begin( );
         itr != differenceStatePerturbedCase.end( ); itr++ ){

        std::cout << "Leg " << itr->first << "\n\n";
        std::cout << "Departure: " << itr->second.first << "\n\n";
        std::cout << "Arrival: " << itr->second.second << "\n\n";

    }



    /// Outputs

    for( std::map< int, std::map< double, Eigen::Vector6d > >::iterator itr = patchedConicsTrajectory.begin( );
         itr != patchedConicsTrajectory.end( ); itr++ ){

        input_output::writeDataMapToTextFile( fullProblemTrajectory[itr->first],
                                              "fullProblemInterplanetaryTrajectory_0_leg_" + std::to_string(itr->first) + ".dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile( fullProblemTrajectoryPerturbedCase[itr->first],
                                              "fullProblemInterplanetaryTrajectory_1_leg_" + std::to_string(itr->first) + ".dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile( patchedConicsTrajectory[itr->first],
                                              "patchedConicsInterplanetaryTrajectory_0_leg_" + std::to_string(itr->first) + ".dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        input_output::writeDataMapToTextFile( patchedConicsTrajectoryPerturbedCase[itr->first],
                                              "patchedConicsInterplanetaryTrajectory_1_leg_" + std::to_string(itr->first) + ".dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

    }


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
