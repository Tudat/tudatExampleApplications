/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      101111    K. Kumar          File created.
 *      110113    K. Kumar          Scenario updated to use latest version of code; added file
 *                                  header and footer.
 *      110202    K. Kumar          Scenario updated to use latest version of code.
 *      110216    K. Kumar          Migrated to applications namespace.
 *      110217    K. Kumar          Function name changed.
 *      110815    K. Kumar          Updated with mass of Asterix.
 *      111024    K. Kumar          Modified to be executable program with main-function as
 *                                  suggested by M. Persson.
 *      120221    K. Kumar          Rewrote application from scratch; now propagates two
 *                                  satellites.
 *
 *    References
 *
 */

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include <boost/bind.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

#include <Tudat/Astrodynamics/Bodies/vehicle.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/forceModel.h>
#include <Tudat/Astrodynamics/Gravitation/gravitationalForceModel.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityField.h>

#include "stateAssembly.h"
#include "stateDerivativeModel.h"

//! Execute example of an Earth-orbiting satellite.
int main( )
{
    using tudat::orbital_element_conversions::semiMajorAxisIndex;
    using tudat::orbital_element_conversions::eccentricityIndex;
    using tudat::orbital_element_conversions::inclinationIndex;
    using tudat::orbital_element_conversions::argumentOfPeriapsisIndex;
    using tudat::orbital_element_conversions::longitudeOfAscendingNodeIndex;
    using tudat::orbital_element_conversions::trueAnomalyIndex;

    using tudat::orbital_element_conversions::xPositionIndex;
    using tudat::orbital_element_conversions::yPositionIndex;
    using tudat::orbital_element_conversions::zPositionIndex;
    using tudat::orbital_element_conversions::xVelocityIndex;
    using tudat::orbital_element_conversions::yVelocityIndex;
    using tudat::orbital_element_conversions::zVelocityIndex;

    using tudat::unit_conversions::convertDegreesToRadians;

    ///////////////////////////////////////////////////////////////////////////

    // Input deck.

    // Set output directory.
    std::string outputDirectory = "";

    // Set simulation start epoch.
    double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    double fixedStepSize = 60.0;

    // Set initial conditions for satellites that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to
    // Cartesian elements.

    // Set Keplerian elements for Asterix.
    Eigen::VectorXd asterixInitialStateInKeplerianElements( 6 );
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0e3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex )
            = convertDegreesToRadians( 139.87 );

    // Set Keplerian elements for Obelix.
    Eigen::VectorXd obelixInitialStateInKeplerianElements( 6 );
    obelixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 12040.6e3;
    obelixInitialStateInKeplerianElements( eccentricityIndex ) = 0.4;
    obelixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( -23.5 );
    obelixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 10.6 );
    obelixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 367.9 );
    obelixInitialStateInKeplerianElements( trueAnomalyIndex )
            = convertDegreesToRadians( 93.4 );

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

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    Eigen::VectorXd asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthCentralGravityField.getGravitationalParameter( ) );

    // Convert Obelix state from Keplerian elements to Cartesian elements.
    Eigen::VectorXd obelixInitialState = convertKeplerianToCartesianElements(
                obelixInitialStateInKeplerianElements,
                earthCentralGravityField.getGravitationalParameter( ) );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Create satellite vehicles and set their respective masses [kg].

    // Declare a new vehicle object for Asterix and set its mass.
    tudat::Vehicle asterix;
    asterix.setMass( 1.0 );

    // Declare a new vehicle object for Obelix and set its mass.
    tudat::Vehicle obelix;
    obelix.setMass( 1.0 );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Create Earth gravitational force models for both satellites.

    // Declare Earth gravitational force models.
    tudat::GravitationalForceModel earthGravitationalForceModelForAsterix(
                &asterix, &earthCentralGravityField );
    tudat::GravitationalForceModel earthGravitationalForceModelForObelix(
                &obelix, &earthCentralGravityField );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Setup forces acting on satellites.

    // Add Earth gravitational forces for Asterix and Obelix.
    using earth_orbiting_satellite_example::StateDerivativeModel;
    StateDerivativeModel::ListOfForces listOfForces;

    std::vector< tudat::ForceModel* > listOfForcesActingOnAsterix;
    listOfForcesActingOnAsterix.push_back( &earthGravitationalForceModelForAsterix );

    listOfForces[ &asterix ] = listOfForcesActingOnAsterix;

    std::vector< tudat::ForceModel* > listOfForcesActingOnObelix;
    listOfForcesActingOnObelix.push_back( &earthGravitationalForceModelForObelix );

    listOfForces[ &obelix ] = listOfForcesActingOnObelix;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Setup initial states for propagation of satellites.

    // Add initial states.
    earth_orbiting_satellite_example::ListOfStates listOfInitialStates;
    listOfInitialStates[ &asterix ] = asterixInitialState;
    listOfInitialStates[ &obelix ] = obelixInitialState;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Setup state derivative model and numerical integrator.

    // Declare state derivative model.
    StateDerivativeModel stateDerivativeModelForEarthSatellites( listOfForces );

    // Declare Runge-Kutta 4 integrator.
    // Since the state derivative function is a member-function, it must be passed by using
    // boost::bind. The "_1" and "_2" in the boost::bind call specifies that the argument list
    // for the computeStateDerivative function takes two arguments (t, x).
    earth_orbiting_satellite_example::AssembledStateWithBodyIndices assembledStateAndBodyIndices;
    assembledStateAndBodyIndices = earth_orbiting_satellite_example::assembleState(
                listOfInitialStates );

    tudat::mathematics::numerical_integrators::RungeKutta4Integrator<
            double, Eigen::VectorXd, Eigen::VectorXd > rungeKutta4(
                boost::bind( &StateDerivativeModel::computeStateDerivative,
                             &stateDerivativeModelForEarthSatellites, _1, _2 ),
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
    std::ofstream asterixPropagationHistoryFile( ( outputDirectory
                                                   + "asterixPropagationHistory.dat" ).c_str( ) );

    std::ofstream obelixPropagationHistoryFile( ( outputDirectory
                                                   + "obelixPropagationHistory.dat" ).c_str( ) );

    // Store state history in propagation history files.
    for ( PropagationHistory::iterator iteratorPropagationHistory = propagationHistory.begin( );
          iteratorPropagationHistory != propagationHistory.end( ); iteratorPropagationHistory++ )
    {
        asterixPropagationHistoryFile
                << iteratorPropagationHistory->first << ", "
                << iteratorPropagationHistory->second[ &asterix ]( xPositionIndex ) << ", "
                << iteratorPropagationHistory->second[ &asterix ]( yPositionIndex ) << ", "
                << iteratorPropagationHistory->second[ &asterix ]( zPositionIndex ) << ", "
                << iteratorPropagationHistory->second[ &asterix ]( xVelocityIndex ) << ", "
                << iteratorPropagationHistory->second[ &asterix ]( yVelocityIndex ) << ", "
                << iteratorPropagationHistory->second[ &asterix ]( zVelocityIndex ) << std::endl;

        obelixPropagationHistoryFile
                << iteratorPropagationHistory->first << ", "
                << iteratorPropagationHistory->second[ &obelix ]( xPositionIndex ) << ", "
                << iteratorPropagationHistory->second[ &obelix ]( yPositionIndex ) << ", "
                << iteratorPropagationHistory->second[ &obelix ]( zPositionIndex ) << ", "
                << iteratorPropagationHistory->second[ &obelix ]( xVelocityIndex ) << ", "
                << iteratorPropagationHistory->second[ &obelix ]( yVelocityIndex ) << ", "
                << iteratorPropagationHistory->second[ &obelix ]( zVelocityIndex ) << std::endl;
    }

    // Close simulation output files.
    asterixPropagationHistoryFile.close( );
    obelixPropagationHistoryFile.close( );

    ///////////////////////////////////////////////////////////////////////////

    return 0;
}
