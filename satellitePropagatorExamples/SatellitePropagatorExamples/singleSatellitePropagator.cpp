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
 *      130916    E.D. Brandon      File created. (used parts of asterixAndObelixPropagator.cpp)
 *
 *    References
 *
 *    Notes
 *      This file serves as a tutorial on how to use different elements of the Tudat Library for
 *      your own simulation. This is the beginner tutorial on numerical propagation of a satellite
 *      orbit. It includes the following important elements that are in the Tudat library:
 *      - the transformation from Keplerian to Cartesian elements.
 *      - the acceleration model.
 *      - the state derivative model.
 *      - the Runge-Kutta 4th-order fixed step-size numerical-integrator.
 *
 */

// ------------------------------------------------------------------------------------------------
// INCLUDE STATEMENTS
// ------------------------------------------------------------------------------------------------

// C++ Standard library
#include <iostream>

// External libraries: Boost
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

// External libraries: Eigen
#include <Eigen/Core>

// Tudat Core library
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

// Tudat library
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

// The simulation's data repository: the body class
#include "body.h"

// ------------------------------------------------------------------------------------------------
// THE MAIN FUNCTION
// ------------------------------------------------------------------------------------------------

//! Execute propagation of orbit of Asterix around the Earth.
int main()
{
    // --------------------------------------------------------------------------------------------
    // USING STATEMENTS (IN ALPHABETICAL ORDER)
    // --------------------------------------------------------------------------------------------

    // The body class.
    using satellite_propagator_examples::Body;

    // Type definition of a shared-pointer to the body class.
    using satellite_propagator_examples::BodyPointer;

    // Keplerian state vector element indices.
    using tudat::basic_astrodynamics::semiMajorAxisIndex;
    using tudat::basic_astrodynamics::eccentricityIndex;
    using tudat::basic_astrodynamics::inclinationIndex;
    using tudat::basic_astrodynamics::argumentOfPeriapsisIndex;
    using tudat::basic_astrodynamics::longitudeOfAscendingNodeIndex;
    using tudat::basic_astrodynamics::trueAnomalyIndex;

    // State vector typedef.
    using tudat::basic_mathematics::Vector6d;

    // Central gravitational acceleration model.
    using tudat::gravitation::CentralGravitationalAccelerationModel3d;
    using tudat::gravitation::CentralGravitationalAccelerationModel3dPointer;

    // Runge-Kutta 4 integrator.
    using tudat::numerical_integrators::RungeKutta4IntegratorXd;

    // Convert Keplerian elements to Cartesian elements.
    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;

    // Cartesian state derivative model.
    using tudat::state_derivative_models::CartesianStateDerivativeModel6d;
    using tudat::state_derivative_models::CartesianStateDerivativeModel6dPointer;

    // Convert degrees to radians.
    using tudat::unit_conversions::convertDegreesToRadians;

    // --------------------------------------------------------------------------------------------
    // INPUT DECK: Constants
    // --------------------------------------------------------------------------------------------

    // Set Earth gravitational parameter [m^3 s^-2].
    const double earthGravitationalParameter = 3.986004415E14;

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 60.0;

    // --------------------------------------------------------------------------------------------
    // INPUT DECK: Initial conditions
    // --------------------------------------------------------------------------------------------

    // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.

    // Set Keplerian elements for Asterix.
    Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // --------------------------------------------------------------------------------------------
    // CONVERT INITIAL STATE FROM KEPLERIAN TO CARTESIAN ELEMENTS
    // --------------------------------------------------------------------------------------------

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    const Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );

    // --------------------------------------------------------------------------------------------
    // CREATE ASTERIX SATELLITE, ACCELERATION MODEL AND STATE DERIVATIVE MODEL
    // --------------------------------------------------------------------------------------------

    // Create a 'Body' data repository for Asterix.
    const BodyPointer asterix = boost::make_shared< Body >( asterixInitialState );

    // Create the gravitational acceleration model for Asterix.
    const CentralGravitationalAccelerationModel3dPointer asterixGravityModel
            = boost::make_shared< CentralGravitationalAccelerationModel3d >(
                boost::bind( &Body::getCurrentPosition, asterix ),
                earthGravitationalParameter );

    // Before creating the Cartesian state derivative model, assign the acceleration model to the
    // AccelerationModelPointerVector, which is an input to the state derivative model.
    CartesianStateDerivativeModel6d::AccelerationModelPointerVector asterixGravity;
    asterixGravity.push_back( asterixGravityModel );

    // Create the Cartesian state derivative model for Asterix.
    const CartesianStateDerivativeModel6dPointer asterixStateDerivativeModel
            = boost::make_shared< CartesianStateDerivativeModel6d >(
                asterixGravity, boost::bind( &Body::setCurrentTimeAndState, asterix, _1, _2 ) );

    // --------------------------------------------------------------------------------------------
    // NUMERICAL INTEGRATION
    // --------------------------------------------------------------------------------------------

    // Set up the numerical integrator: the Runge-Kutta 4 integrator.
    // Since the state derivative function is a member-function, it must be passed by using
    // boost::bind. The "_1" and "_2" in the boost::bind call specifies that the argument list
    // for the computeStateDerivative function takes two arguments (t, x).
    RungeKutta4IntegratorXd rungeKutta4(
                boost::bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
                             asterixStateDerivativeModel, _1, _2 ),
                0.0, asterixInitialState );

    // Execute the simulation and compute the final state.
    Vector6d finalIntegratedState = rungeKutta4.integrateTo( simulationEndEpoch, fixedStepSize );

    // Print the position (in km) and the velocity (in km/s) at t = 0.
    std::cout << "Single Earth-Orbiting Satellite Example." << std::endl <<
                 "The initial position vector of Asterix is [km]:" << std::endl <<
                 asterixInitialState.segment( 0, 3 ) / 1E3 << std::endl <<
                 "The initial velocity vector of Asterix is [km/s]:" << std::endl <<
                 asterixInitialState.segment( 3, 3 ) / 1E3 << std::endl;

    // Print the position (in km) and the velocity (in km/s) at t = 86400.
    std::cout << "After " << simulationEndEpoch <<
                 " seconds, the position vector of Asterix is [km]:" << std::endl <<
                 finalIntegratedState.segment( 0, 3 ) / 1E3 << std::endl <<
                 "And the velocity vector of Asterix is [km/s]:" << std::endl <<
                 finalIntegratedState.segment( 3, 3 ) / 1E3 << std::endl;

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
