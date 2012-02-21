/*!   Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120221    K. Kumar          File created.
 *
 *    References
 *
 */

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

#include <Tudat/Astrodynamics/Bodies/vehicle.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityField.h>
#include <Tudat/Astrodynamics/Gravitation/gravitationalForceModel.h>

#include "stateAssembly.h"
#include "stateDerivativeModel.h"

//! Execute simulation of Galileo constellation around the Earth.
int main( )
{
    using tudat::unit_conversions::convertDegreesToRadians;

    using tudat::orbital_element_conversions::xPositionIndex;
    using tudat::orbital_element_conversions::yPositionIndex;
    using tudat::orbital_element_conversions::zPositionIndex;
    using tudat::orbital_element_conversions::xVelocityIndex;
    using tudat::orbital_element_conversions::yVelocityIndex;
    using tudat::orbital_element_conversions::zVelocityIndex;

    ///////////////////////////////////////////////////////////////////////////

    // Input deck.

    // Set output directory.
    std::string outputDirectory = "";

    // Set simulation start epoch.
    double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    double simulationEndEpoch = 0.1 * tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    double fixedStepSize = 30.0;

    // Set number of satellites in constellation.
    const unsigned int NUMBER_OF_SATELLITES = 30;

    // Set orbital parameters of Galileo constellation.
    double semiMajorAxis = 21000.0e3;                                           // [km]
    double eccentricity = 0.0;                                                  // [-]
    double inclination = convertDegreesToRadians( 60.0 );                       // [rad]
    double argumentOfPeriapsis = 0.0;                                           // [rad]
    double longitudeOfAscendingNodeSpacing = convertDegreesToRadians( 120.0 );  // [rad]
    double trueAnomalySpacing = convertDegreesToRadians( 36.0 );                // [rad]

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set initial conditions for Galileo satellites that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to
    // Cartesian elements. They are stored in a matrix.

    // Declare size of state.
    const unsigned int SIZE_OF_STATE = 6;

    // Set Keplerian elements for Galileo satellites.
    Eigen::MatrixXd initialConditionsInKeplerianElements( NUMBER_OF_SATELLITES, SIZE_OF_STATE );

    // Set semiMajorAxis.
    initialConditionsInKeplerianElements.col( 0 )
            = Eigen::MatrixXd::Constant( NUMBER_OF_SATELLITES, 1, semiMajorAxis );

    // Set eccentricity.
    initialConditionsInKeplerianElements.col( 1 )
            = Eigen::MatrixXd::Constant( NUMBER_OF_SATELLITES, 1, eccentricity );

    // Set inclination.
    initialConditionsInKeplerianElements.col( 2 )
            = Eigen::MatrixXd::Constant( NUMBER_OF_SATELLITES, 1, inclination );

    // Set argument of periapsis.
    initialConditionsInKeplerianElements.col( 3 )
            = Eigen::MatrixXd::Constant( NUMBER_OF_SATELLITES, 1, argumentOfPeriapsis );

    // Set longitude of ascending node.
    initialConditionsInKeplerianElements.block( 0, 4, 10, 1 )
            = Eigen::MatrixXd::Constant( 10, 1, 0.0 );
    initialConditionsInKeplerianElements.block( 10, 4, 10, 1 )
            = Eigen::MatrixXd::Constant( 10, 1, 1.0 * longitudeOfAscendingNodeSpacing );
    initialConditionsInKeplerianElements.block( 20, 4, 10, 1 )
            = Eigen::MatrixXd::Constant( 10, 1, 2.0 * longitudeOfAscendingNodeSpacing );

    // Set true anomaly.
    Eigen::VectorXd trueAnomalySpacingIntegers( 10 );
    trueAnomalySpacingIntegers << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    initialConditionsInKeplerianElements.block( 0, 5, 10, 1 )
            = Eigen::MatrixXd::Constant( 10, 1, trueAnomalySpacing ).array( )
            * trueAnomalySpacingIntegers.array( );
    initialConditionsInKeplerianElements.block( 10, 5, 10, 1 )
            = Eigen::MatrixXd::Constant( 10, 1, trueAnomalySpacing ).array( )
            * trueAnomalySpacingIntegers.array( );
    initialConditionsInKeplerianElements.block( 20, 5, 10, 1 )
            = Eigen::MatrixXd::Constant( 10, 1, trueAnomalySpacing ).array( )
            * trueAnomalySpacingIntegers.array( );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Declare Earth environment.

    // Declare pre-defined Earth central gravity field.
    tudat::CentralGravityField earthCentralGravityField;
    earthCentralGravityField.setPredefinedCentralGravityFieldSettings(
                tudat::CentralGravityField::earth );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Convert initial states from Keplerian to Cartesian elements.

    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;

    // Set Cartesian elements for Galileo satellites.
    Eigen::MatrixXd initialConditions( NUMBER_OF_SATELLITES, SIZE_OF_STATE );

    for ( unsigned int i = 0; i < initialConditionsInKeplerianElements.rows( ); i++ )
    {
        // Convert state from Keplerian elements to Cartesian elements.
        initialConditions.row( i ) = convertKeplerianToCartesianElements(
                    initialConditionsInKeplerianElements.row( i ),
                    earthCentralGravityField.getGravitationalParameter( ) );
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Create satellite vehicles and set their respective masses [kg].

    std::vector< tudat::Vehicle > vehicles( NUMBER_OF_SATELLITES );

    for ( unsigned int i = 0; i < NUMBER_OF_SATELLITES; i++ )
    {
        vehicles.at( i ).setMass( 1.0 );
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Create Earth gravitational force models for satellites.

    // Declare Earth gravitational force models.
    std::vector< tudat::GravitationalForceModel > earthGravitationalForceModels;

    for ( unsigned int i = 0; i < NUMBER_OF_SATELLITES; i++ )
    {
        tudat::GravitationalForceModel earthGravitationalForceModel(
                    &vehicles.at( i ), &earthCentralGravityField );
        earthGravitationalForceModels.push_back( earthGravitationalForceModel );
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Setup forces acting on satellites.

    // Add Earth gravitational forces for Asterix and Obelix.
    using earth_orbiting_satellite_example::StateDerivativeModel;
    StateDerivativeModel::ListOfForces listOfForces;

    for ( unsigned int i = 0; i < NUMBER_OF_SATELLITES; i++ )
    {
        std::vector< tudat::ForceModel* > forceModels;
        forceModels.push_back( &earthGravitationalForceModels.at( i ) );
        listOfForces[ &vehicles.at( i ) ] = forceModels;
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Setup initial states for propagation of satellites.

    // Add initial states.
    earth_orbiting_satellite_example::ListOfStates listOfInitialStates;

    for ( unsigned int i = 0; i < NUMBER_OF_SATELLITES; i++ )
    {
        listOfInitialStates[ &vehicles.at( i ) ] = initialConditions.row( i );
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Setup state derivative model and numerical integrator.

    // Declare state derivative model.
    StateDerivativeModel stateDerivativeModelForSatellites( listOfForces );

    // Declare Runge-Kutta 4 integrator.
    // Since the state derivative function is a member-function, it must be passed by using
    // boost::bind. The "_1" and "_2" in the boost::bind call specifies that the argument list
    // for the computeStateDerivative function takes two arguments (t, x).
    earth_orbiting_satellite_example::AssembledStateWithBodyIndices assembledStateAndBodyIndices;
    assembledStateAndBodyIndices = earth_orbiting_satellite_example::assembleState(
                listOfInitialStates );

    tudat::mathematics::numerical_integrators::RungeKutta4IntegratorXd rungeKutta4(
                boost::bind( &StateDerivativeModel::computeStateDerivative,
                             &stateDerivativeModelForSatellites, _1, _2 ),
                0.0, assembledStateAndBodyIndices.first );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Execute simulation.

    // Set running time, updated after each step that the numerical integrator takes.
    double runningTime = simulationStartEpoch;

    // Declare propagation history to store state history of satellites.
    using earth_orbiting_satellite_example::ListOfStates;
    typedef std::map< double, ListOfStates > PropagationHistory;
    PropagationHistory propagationHistory;

    // Set initial states in propagation history.
    propagationHistory[ 0.0 ] = listOfInitialStates;

    // Execute simulation from start to end epoch and save intermediate states in propagation
    // history.
    while ( runningTime < simulationEndEpoch )
    {
        // Execute integration step and store state at end.
        Eigen::VectorXd integratedState = rungeKutta4.performIntegrationStep( fixedStepSize );

        // Update running time to value of current time.
        runningTime = rungeKutta4.getCurrentIndependentVariable( );

        // Disassemble state into states per body, and store in propagation history.
        propagationHistory[ runningTime ]
                = earth_orbiting_satellite_example::disassembleState(
                    std::make_pair( integratedState, assembledStateAndBodyIndices.second ) );
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Write results to file.

    // Declare results file for propagation history for satellites.
    std::vector< boost::shared_ptr< std::ofstream > > propagationHistoryDataFiles;

    for ( unsigned int i = 0; i < NUMBER_OF_SATELLITES; i++ )
    {
        std::stringstream propagationHistoryDataFilename;
        propagationHistoryDataFilename << outputDirectory << "galileoSatellite" << i + 1 << ".dat";
        propagationHistoryDataFiles.push_back(
                    boost::shared_ptr< std::ofstream >(
                        new std::ofstream( propagationHistoryDataFilename.str( ).c_str( ) ) ) );
    }

    // Store state history in propagation history files.
    for ( PropagationHistory::iterator iteratorPropagationHistory = propagationHistory.begin( );
          iteratorPropagationHistory != propagationHistory.end( ); iteratorPropagationHistory++ )
    {
        for ( unsigned int i = 0; i < NUMBER_OF_SATELLITES; i++ )
        {
            *propagationHistoryDataFiles.at( i )
                    << iteratorPropagationHistory->first << ", "
                    << iteratorPropagationHistory->second[ &vehicles.at( i ) ]( xPositionIndex ) << ", "
                    << iteratorPropagationHistory->second[ &vehicles.at( i ) ]( yPositionIndex ) << ", "
                    << iteratorPropagationHistory->second[ &vehicles.at( i ) ]( zPositionIndex ) << ", "
                    << iteratorPropagationHistory->second[ &vehicles.at( i ) ]( xVelocityIndex ) << ", "
                    << iteratorPropagationHistory->second[ &vehicles.at( i ) ]( yVelocityIndex ) << ", "
                    << iteratorPropagationHistory->second[ &vehicles.at( i ) ]( zVelocityIndex ) << std::endl;
        }
    }

    // Close simulation output files.
//    propagationHistoryDataFile.close( );

    ///////////////////////////////////////////////////////////////////////////

    return 0;

}
