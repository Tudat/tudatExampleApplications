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
 *      120502    K. Kumar          Updated code to use shared pointers.
 *      121030    K. Kumar          Updated code to use new state derivative models.
 *
 *    References
 *
 *    Notes
 *
 */

#include <fstream>
#include <limits>
#include <string>
#include <utility>

#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Gravitation/centralJ2J3J4GravityModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h>

#include <Tudat/InputOutput/basicInputOutput.h>

#include "body.h"

//! Execute example of an Earth-orbiting satellite.
int main( )
{
    using namespace satellite_example;

    using tudat::basic_astrodynamics::AccelerationModel3dPointer;
    using tudat::basic_astrodynamics::semiMajorAxisIndex;
    using tudat::basic_astrodynamics::eccentricityIndex;
    using tudat::basic_astrodynamics::inclinationIndex;
    using tudat::basic_astrodynamics::argumentOfPeriapsisIndex;
    using tudat::basic_astrodynamics::longitudeOfAscendingNodeIndex;
    using tudat::basic_astrodynamics::trueAnomalyIndex;
    using tudat::basic_astrodynamics::xCartesianPositionIndex;
    using tudat::basic_astrodynamics::yCartesianPositionIndex;
    using tudat::basic_astrodynamics::zCartesianPositionIndex;
    using tudat::basic_astrodynamics::xCartesianVelocityIndex;
    using tudat::basic_astrodynamics::yCartesianVelocityIndex;
    using tudat::basic_astrodynamics::zCartesianVelocityIndex;

    using tudat::gravitation::CentralJ2J3J4GravitationalAccelerationModel3d;
    using tudat::gravitation::CentralJ2J3J4GravitationalAccelerationModel3dPointer;

    using tudat::input_output::writeDataMapToTextFile;

    using tudat::mathematics::numerical_integrators::RungeKutta4Integrator;

    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;

    using tudat::state_derivative_models::CartesianStateDerivativeModel6d;
    using tudat::state_derivative_models::CartesianStateDerivativeModel6dPointer;
    using tudat::state_derivative_models::CompositeStateDerivativeModel;

    using tudat::unit_conversions::convertDegreesToRadians;

    typedef Eigen::Matrix< double, 12, 1 > Vector12d;
    typedef CompositeStateDerivativeModel< double, Vector12d, Vector6d >
            CompositionStateDerivativeModel12d;

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
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Set Keplerian elements for Obelix.
    Eigen::VectorXd obelixInitialStateInKeplerianElements( 6 );
    obelixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 12040.6e3;
    obelixInitialStateInKeplerianElements( eccentricityIndex ) = 0.4;
    obelixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( -23.5 );
    obelixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 10.6 );
    obelixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 367.9 );
    obelixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 93.4 );

    // Set Earth gravitational parameter [m^3 s^-2].
    const double earthGravitationalParameter = 3.986004415e14;

    // Set spherical harmonics zonal term coefficients.
    const double earthJ2 = 0.0010826269;
    const double earthJ3 = -0.0000025323;
    const double earthJ4 = -0.0000016204;

    // Set equatorial radius of Earth [m].
    const double earthEquatorialRadius = 6378.1363e3;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Convert initial states from Keplerian to Cartesian elements.

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    Eigen::VectorXd asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );

    // Convert Obelix state from Keplerian elements to Cartesian elements.
    Eigen::VectorXd obelixInitialState = convertKeplerianToCartesianElements(
                obelixInitialStateInKeplerianElements,
                earthGravitationalParameter );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Create Asterix and Obelix satellites, and gravitational acceleration models.

    // Create Asterix and set initial state and epoch.
    BodyPointer asterix = boost::make_shared< Body >(
                asterixInitialStateInKeplerianElements, 0.0 );

    // Create gravitational acceleration model for asterix.
    CentralJ2J3J4GravitationalAccelerationModel3dPointer asterixGravityModel
            = boost::make_shared< CentralJ2J3J4GravitationalAccelerationModel3d >(
                boost::bind( &Body::getCurrentPosition, asterix ),
                earthGravitationalParameter, earthEquatorialRadius,
                earthJ2, earthJ3, earthJ4 );

    // Create Cartesian state derivative model for asterix.
    CartesianStateDerivativeModel6d::AccelerationModelPointerVector asterixGravity
            = boost::assign::list_of( asterixGravityModel );

    // Create Obelix and set initial state and epoch.
    BodyPointer obelix = boost::make_shared< Body >(
                asterixInitialStateInKeplerianElements, 0.0 );

    // Create gravitational acceleration model for obelix.
    CentralJ2J3J4GravitationalAccelerationModel3dPointer obelixGravityModel
            = boost::make_shared< CentralJ2J3J4GravitationalAccelerationModel3d >(
                boost::bind( &Body::getCurrentPosition, asterix ),
                earthGravitationalParameter, earthEquatorialRadius,
                earthJ2, earthJ3, earthJ4 );

    // Create Cartesian state derivative model model for obelix.
    CartesianStateDerivativeModel6d::AccelerationModelPointerVector obelixGravity
            = boost::assign::list_of( obelixGravityModel );

    // Add Asterix and Obelix to list of satellites.
    ListOfSatellites satellites;
    satellites[ asterix ] = asterixGravity;
    satellites[ obelix ] = obelixGravity;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Create Cartesian state derivative models for Asterix and Obelix.

    // Create Cartesian state derivative model for Asterix.
    CartesianStateDerivativeModel6dPointer asterixStateDerivativeModel
            = boost::make_shared< CartesianStateDerivativeModel6d >(
                asterixGravity, boost::bind( &Body::setCurrentTimeAndState, asterix, _1, _2 ) );

    // Create Cartesian state derivative model for Obelix.
    CartesianStateDerivativeModel6dPointer obelixStateDerivativeModel
            = boost::make_shared< CartesianStateDerivativeModel6d >(
                obelixGravity, boost::bind( &Body::setCurrentTimeAndState, obelix, _1, _2 ) );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Construct composite state derivative model for Asterix and Obelix.

    // Create state derivative model map and bind Asterix and Obelix, in that order.
    CompositionStateDerivativeModel12d::VectorStateDerivativeModelMap stateDerivativeModelMap;

    stateDerivativeModelMap[ std::make_pair( 0, 6 ) ]
            = boost::bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
                           asterixStateDerivativeModel, _1, _2 );

    stateDerivativeModelMap[ std::make_pair( 6, 6 ) ]
            = boost::bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
                           asterixStateDerivativeModel, _1, _2 );

    // Create data updater.
    DataUpdater updater( satellites );

    // Create composite state derivative model.
    boost::shared_ptr< CompositionStateDerivativeModel12d > stateDerivativeModel
            = boost::make_shared< CompositionStateDerivativeModel12d >(
                stateDerivativeModelMap,
                boost::bind( &DataUpdater::updateBodyData, updater, _1, _2 ) );

    ///////////////////////////////////////////////////////////////////////////

    // Set up numerical integrator and execute simulation.

    // Declare Runge-Kutta 4 integrator.
    // Since the state derivative function is a member-function, it must be passed by using
    // boost::bind. The "_1" and "_2" in the boost::bind call specifies that the argument list
    // for the computeStateDerivative function takes two arguments (t, x).
    RungeKutta4Integrator< double, Vector12d, Vector12d > rungeKutta4(
                boost::bind( &CompositionStateDerivativeModel12d::computeStateDerivative,
                             stateDerivativeModel, _1, _2 ),
                0.0, ( Eigen::VectorXd( 12 ) << asterixInitialState,
                       obelixInitialState ).finished( ) );

    // Set running time, updated after each step that the numerical integrator takes.
    double runningTime = simulationStartEpoch;

    // Declare propagation history to store state history of satellites.
    PropagationHistory asterixPropagationHistory;
    PropagationHistory obelixPropagationHistory;

    // Set initial states in propagation history.
    asterixPropagationHistory[ 0.0 ] = asterixInitialState;
    obelixPropagationHistory[ 0.0 ] = obelixInitialState;

    // Execute simulation from start to end epoch and save intermediate states in propagation
    // history.
    while ( runningTime < simulationEndEpoch )
    {
        // Execute integration step and store state at end.
        Vector12d integratedState = rungeKutta4.performIntegrationStep( fixedStepSize );

        // Update running time to value of current time.
        runningTime = rungeKutta4.getCurrentIndependentVariable( );

        // Disassemble state into states per body, and store in propagation history.
        asterixPropagationHistory[ runningTime ] = integratedState.segment( 0, 6 );
        obelixPropagationHistory[ runningTime ] = integratedState.segment( 6, 6 );
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Write results to file.

    // Write Asterix propagation history to file.
    writeDataMapToTextFile( asterixPropagationHistory,
                            "asterixPropagationHistory.dat",
                            outputDirectory,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    // Write obelix propagation history to file.
    writeDataMapToTextFile( obelixPropagationHistory,
                            "obelixPropagationHistory.dat",
                            outputDirectory,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    ///////////////////////////////////////////////////////////////////////////

    return 0;
}
