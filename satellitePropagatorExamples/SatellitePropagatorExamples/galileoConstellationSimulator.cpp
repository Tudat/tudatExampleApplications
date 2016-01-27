/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120221    K. Kumar          File created.
 *      120502    K. Kumar          Updated code to use shared pointers.
 *      121030    K. Kumar          Updated code to use new state derivative models.
 *      130107    S. Billemont      Fixed bugs in set up of initial conditions.
 *      130107    K. Kumar          Updated license in file header.
 *      130225    K. Kumar          Updated gravitational acceleration model references; renamed
 *                                  file; made variables const-correct.
 *
 *    References
 *
 *    Notes
 *
 */

#include <boost/assign/list_of.hpp>
#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/multi_array.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Gravitation/centralJ2J3J4GravityModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "SatellitePropagatorExamples/body.h"

//! Execute simulation of Galileo constellation around the Earth.
int main( )
{
    using namespace satellite_propagator_examples;

    using tudat::orbital_element_conversions::xCartesianPositionIndex;
    using tudat::orbital_element_conversions::yCartesianPositionIndex;
    using tudat::orbital_element_conversions::zCartesianPositionIndex;
    using tudat::orbital_element_conversions::xCartesianVelocityIndex;
    using tudat::orbital_element_conversions::yCartesianVelocityIndex;
    using tudat::orbital_element_conversions::zCartesianVelocityIndex;

    using tudat::basic_mathematics::Vector6d;

    using tudat::gravitation::CentralJ2J3J4GravitationalAccelerationModel;

    using tudat::input_output::DoubleKeyTypeVectorXdValueTypeMap;
    using tudat::input_output::writeDataMapToTextFile;

    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;

    using tudat::numerical_integrators::RungeKutta4Integrator;

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
    const std::string outputDirectory = "";
    if( outputDirectory == "" )
    {
        std::cerr<<"Error, output directory not specified (modify outputDirectory variable to "<<
                   " required output directory)."<<std::endl;
    }

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
            = 2.0 * tudat::mathematical_constants::PI / numberOfPlanes;                          // [rad]
    const double trueAnomalySpacing
            = 2.0 * tudat::mathematical_constants::PI / numberOfSatellitesPerPlane;              // [rad]

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
    for ( unsigned int i = 0; i < numberOfPlanes; i++ )

    {
        initialConditionsInKeplerianElements.block( 4, i * numberOfSatellitesPerPlane,
                                                    1, numberOfSatellitesPerPlane )
                = Eigen::MatrixXd::Constant( 1, numberOfSatellitesPerPlane,
                                             i * longitudeOfAscendingNodeSpacing );
    }

    // Set true anomaly.
    Eigen::RowVectorXd trueAnomalySpacingIntegers( numberOfSatellitesPerPlane );

    for ( unsigned int i = 0; i < numberOfSatellitesPerPlane; i++ )
    {
        trueAnomalySpacingIntegers( i ) =  i * 1.0;
    }

    for ( unsigned int i = 0; i < numberOfPlanes; i++ )
    {
        initialConditionsInKeplerianElements.block( 5, i * numberOfSatellitesPerPlane,
                                                    1, numberOfSatellitesPerPlane )
                = Eigen::MatrixXd::Constant( 1, numberOfSatellitesPerPlane,
                                             trueAnomalySpacing ).array( )
                * trueAnomalySpacingIntegers.array( );
    }

    // Convert initial conditions to Cartesian elements.
    Eigen::MatrixXd initialConditions( sizeOfState, numberOfSatellites );

    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
		Vector6d initKepl = initialConditionsInKeplerianElements.col( i ).cast< double >();
        initialConditions.col( i ) = convertKeplerianToCartesianElements(
                    initKepl, static_cast< double >(earthGravitationalParameter) );
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
		
        const CartesianStateDerivativeModel6d::AccelerationModelPointerVector gravityModel = 
		    boost::assign::list_of( boost::make_shared< CentralJ2J3J4GravitationalAccelerationModel >(
                        boost::bind( &Body::getCurrentPosition, satellite ),
                        earthGravitationalParameter, earthEquatorialRadius,
                        earthJ2, earthJ3, earthJ4 ) );
		
        satellites[ satellite ] = gravityModel;
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
    std::vector< DoubleKeyTypeVectorXdValueTypeMap > allSatellitesPropagationHistory(
                numberOfSatellites );

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
