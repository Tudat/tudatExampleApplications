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
 *      120221    K. Kumar          File created.
 *      120502    K. Kumar          Updated code to use shared pointers.
 *      121030    K. Kumar          Updated code to use new state derivative models.
 *
 *    References
 *
 *    Notes
 *
 */

#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <TudatCore/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Gravitation/centralJ2J3J4GravityModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h>

#include <Tudat/InputOutput/basicInputOutput.h>

#include "body.h"

//! Execute simulation of Galileo constellation around the Earth.
int main( )
{
    using namespace satellite_example;

    using tudat::basic_astrodynamics::xCartesianPositionIndex;
    using tudat::basic_astrodynamics::yCartesianPositionIndex;
    using tudat::basic_astrodynamics::zCartesianPositionIndex;
    using tudat::basic_astrodynamics::xCartesianVelocityIndex;
    using tudat::basic_astrodynamics::yCartesianVelocityIndex;
    using tudat::basic_astrodynamics::zCartesianVelocityIndex;

    using tudat::gravitation::CentralJ2J3J4GravitationalAccelerationModel3d;

    using tudat::input_output::writeDataMapToTextFile;

    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;

    using tudat::mathematics::numerical_integrators::RungeKutta4Integrator;

    using tudat::state_derivative_models::CartesianStateDerivativeModel6d;
    using tudat::state_derivative_models::CartesianStateDerivativeModel6dPointer;

    using tudat::unit_conversions::convertDegreesToRadians;

    typedef std::vector< CartesianStateDerivativeModel6dPointer > ListOfStateDerivativeModels;
    typedef tudat::state_derivative_models::CompositeStateDerivativeModel<
            double, Eigen::VectorXd, Vector6d > CompositeStateDerivativeModelXd;
    typedef boost::shared_ptr< CompositeStateDerivativeModelXd >
            CompositeStateDerivativeModelXdPointer;

    ///////////////////////////////////////////////////////////////////////////

    // Input deck.

    // Set output directory.
    std::string outputDirectory = "";

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 1.0 * tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 30.0;

    // Set number of satellites in constellation.
    const unsigned int numberOfSatellites = 30;

    // Set number of planes in constellation.
    const unsigned int numberOfPlanes = 3;

    // Set number of satellites per plane in constellation.
    const unsigned int numberOfSatellitesPerPlane = numberOfSatellites / numberOfPlanes;

    // Define Earth parameters.
    const double earthGravitationalParameter = 3.986004415e14;
    const double earthJ2 = 0.0010826269;
    const double earthJ3 = -0.0000025323;
    const double earthJ4 = -0.0000016204;
    const double earthEquatorialRadius = 6378.1363e3;

    // Set orbital parameters of Galileo constellation.
    const double semiMajorAxis = 23222.0e3 + 6378.1e3;                                // [km]
    const double eccentricity = 0.0;                                                  // [-]
    const double inclination = convertDegreesToRadians( 56.0 );                       // [rad]
    const double argumentOfPeriapsis = 0.0;                                           // [rad]
    const double longitudeOfAscendingNodeSpacing
            = 2.0 * tudat::mathematics::PI / numberOfPlanes;                          // [rad]
    const double trueAnomalySpacing
            = 2.0 * tudat::mathematics::PI / numberOfSatellitesPerPlane;              // [rad]

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set initial conditions for Galileo satellites that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to
    // Cartesian elements. They are stored in a matrix.

    // Declare size of state.
    const unsigned int sizeOfState = 6;

    // Set Keplerian elements for Galileo satellites.
    Eigen::MatrixXd initialConditionsInKeplerianElements( sizeOfState, numberOfSatellites );

    // Set semiMajorAxis.
    initialConditionsInKeplerianElements.row( 0 )
            = Eigen::MatrixXd::Constant( 1, numberOfSatellites, semiMajorAxis );

    // Set eccentricity.
    initialConditionsInKeplerianElements.row( 1 )
            = Eigen::MatrixXd::Constant( 1, numberOfSatellites, eccentricity );

    // Set inclination.
    initialConditionsInKeplerianElements.row( 2 )
            = Eigen::MatrixXd::Constant( 1, numberOfSatellites, inclination );

    // Set argument of periapsis.
    initialConditionsInKeplerianElements.row( 3 )
            = Eigen::MatrixXd::Constant( 1, numberOfSatellites, argumentOfPeriapsis );

    // Set longitude of ascending node.
    for ( unsigned int i = 0; i <= numberOfPlanes; i++ )

    {
        initialConditionsInKeplerianElements.block( 4, i * numberOfSatellitesPerPlane,
                                                    1, numberOfSatellitesPerPlane )
                = Eigen::VectorXd::Constant( 1, numberOfSatellitesPerPlane,
                                             i * longitudeOfAscendingNodeSpacing );
    }

    // Set true anomaly.
    Eigen::RowVectorXd trueAnomalySpacingIntegers( numberOfSatellitesPerPlane );

    for ( unsigned int i = 0; i < numberOfSatellitesPerPlane; i++ )
    {
        trueAnomalySpacingIntegers( i ) =  i * 1.0;
    }

    for ( unsigned int i = 0; i <= numberOfSatellitesPerPlane; i++ )
    {
        initialConditionsInKeplerianElements.block( 5, i * numberOfSatellitesPerPlane,
                                                    1, numberOfSatellitesPerPlane )
                = Eigen::VectorXd::Constant( 1, numberOfSatellitesPerPlane,
                                             trueAnomalySpacing ).array( )
                * trueAnomalySpacingIntegers.array( );
    }

    // Convert initial conditions to Cartesian elements.
    Eigen::MatrixXd initialConditions( sizeOfState, numberOfSatellites );

    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        initialConditions.col( i ) = convertKeplerianToCartesianElements(
                    initialConditionsInKeplerianElements.col( i ), earthGravitationalParameter );
    }

//    // DEBUG.
//    std::cout << initialConditionsInKeplerianElements << std::endl;
//    std::cout << initialConditions << std::endl;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Create list of satellites and acceleration models per satellite.
    ListOfSatellites satellites;

    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        BodyPointer satellite = boost::make_shared< Body >( initialConditions.col( i ), 0.0 );

        satellites[ satellite ]
                = boost::assign::list_of(
                    boost::make_shared< CentralJ2J3J4GravitationalAccelerationModel3d >(
                        boost::bind( &Body::getCurrentPosition, satellite ),
                        earthGravitationalParameter, earthEquatorialRadius,
                        earthJ2, earthJ3, earthJ4 ) );
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Create list of Cartesian state derivative models.
    ListOfStateDerivativeModels stateDerivativeModels;

    for ( ListOfSatellites::iterator satelliteIterator = satellites.begin( );
          satelliteIterator != satellites.end( ); satelliteIterator++ )
    {
        stateDerivativeModels.push_back(
                    boost::make_shared< CartesianStateDerivativeModel6d >(
                        satelliteIterator->second,
                        boost::bind( &Body::setCurrentTimeAndState,
                                     satelliteIterator->first, _1, _2 ) ) );
    }

    // Create state derivative model map.
    CompositeStateDerivativeModelXd::VectorStateDerivativeModelMap stateDerivativeModelMap;

    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        stateDerivativeModelMap[ std::make_pair( i * sizeOfState, sizeOfState ) ]
                = boost::bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
                               stateDerivativeModels.at( i ), _1, _2 );
    }

    // Create data updater.
    DataUpdater updater( satellites );

    // Create composite state derivative model.
    CompositeStateDerivativeModelXdPointer compositeStateDerivativeModel
            = boost::make_shared< CompositeStateDerivativeModelXd >(
                stateDerivativeModelMap,
                boost::bind( &DataUpdater::updateBodyData, updater, _1, _2 ) );

//    // DEBUG.
//    std::cout << compositeStateDerivativeModel->computeStateDerivative(
//                     0.0, initialConditions ) << std::endl;

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Set up numerical integrator and execute simulation..

    // Integrator initial state.
    Eigen::VectorXd initialState = Eigen::VectorXd( numberOfSatellites * sizeOfState );

    unsigned int counter = 0;

    for ( ListOfSatellites::iterator satelliteIterator = satellites.begin( );
          satelliteIterator != satellites.end( ); satelliteIterator++ )
    {
        initialState.segment( counter * sizeOfState, sizeOfState )
                = satelliteIterator->first->getCurrentState( );
        counter++;
    }

//    // DEBUG.
//    std::cout << initialState << std::endl;

    // Declare Runge-Kutta 4 integrator.
    // Since the state derivative function is a member-function, it must be passed by using
    // boost::bind. The "_1" and "_2" in the boost::bind call specifies that the argument list
    // for the computeStateDerivative function takes two arguments (t, x).
    RungeKutta4Integrator< double, Eigen::VectorXd, Eigen::VectorXd > rungeKutta4(
                boost::bind( &CompositeStateDerivativeModelXd::computeStateDerivative,
                             compositeStateDerivativeModel, _1, _2 ),
                0.0, initialState );

    // Numerically integrate motion of satellites.

    // Set current simulation epoch to start.
    double currentEpoch = simulationStartEpoch;

    // Declare vector of propagation histories for satellites.
    std::vector< PropagationHistory > allSatellitesPropagationHistory( numberOfSatellites );

    // Execute simulation from start to end epoch and save intermediate states in vector of
    // propagation histories.
    while ( currentEpoch < simulationEndEpoch )
    {
        // Execute integration step and store composite state at end.
        Eigen::VectorXd currentCompositeState = rungeKutta4.integrateTo(
                    currentEpoch + fixedStepSize, 1.0 );

        // Update current epoch.
        currentEpoch += fixedStepSize;

        // Disassemble state into states per body, and store in propagation histories.
        for ( unsigned int i = 0; i < numberOfSatellites; i++ )
        {
            allSatellitesPropagationHistory.at( i )[ rungeKutta4.getCurrentIndependentVariable( ) ]
                    = currentCompositeState.segment( i * sizeOfState, sizeOfState );

        }
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    // Write results to file.

    // Loop over all satellites.
    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        // Set filename for output data.
        std::stringstream outputFilename;
        outputFilename << "galileoSatellite" << i + 1 << ".dat";

        // Write propagation history to file.
        writeDataMapToTextFile( allSatellitesPropagationHistory.at( i ),
                                outputFilename.str( ),
                                outputDirectory,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                ", " );
    }

    ///////////////////////////////////////////////////////////////////////////

    return 0;
}
