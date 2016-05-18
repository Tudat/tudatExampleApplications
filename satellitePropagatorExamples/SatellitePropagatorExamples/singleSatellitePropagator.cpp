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


// External libraries: Boost
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

// Tudat Core library
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

// Tudat library
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <Tudat/InputOutput/basicInputOutput.h>

#include <Tudat/Astrodynamics/Propagators/dynamicsSimulator.h>
#include <Tudat/External/SpiceInterface/spiceInterface.h>
#include <Tudat/SimulationSetup/body.h>
#include <Tudat/SimulationSetup/createBodies.h>
#include <Tudat/SimulationSetup/createAccelerationModels.h>

// C++ Standard library
#include <iostream>

// External libraries: Eigen
#include <Eigen/Core>


// ------------------------------------------------------------------------------------------------
// THE MAIN FUNCTION
// ------------------------------------------------------------------------------------------------




//! Execute propagation of orbit of Asterix around the Earth.
int main()
{


    // --------------------------------------------------------------------------------------------
    // USING STATEMENTS
    // --------------------------------------------------------------------------------------------ter;

    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;


    // Keplerian state vector element indices.
    using tudat::orbital_element_conversions::semiMajorAxisIndex;
    using tudat::orbital_element_conversions::eccentricityIndex;
    using tudat::orbital_element_conversions::inclinationIndex;
    using tudat::orbital_element_conversions::argumentOfPeriapsisIndex;
    using tudat::orbital_element_conversions::longitudeOfAscendingNodeIndex;
    using tudat::orbital_element_conversions::trueAnomalyIndex;

    // State vector typedef.
    using tudat::basic_mathematics::Vector6d;

    // Central gravitational acceleration model.
    using tudat::gravitation::CentralGravitationalAccelerationModel3d;
    using tudat::gravitation::CentralGravitationalAccelerationModel3dPointer;

    // Runge-Kutta 4 integrator.
    using tudat::numerical_integrators::RungeKutta4IntegratorXd;

    // Convert Keplerian elements to Cartesian elements.
    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;

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

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Define body settings for simulation.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                basic_mathematics::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralBodyMap;
    Eigen::VectorXd systemInitialState = asterixInitialState;

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationMap[  "Asterix" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );
    centralBodyMap[  "Asterix" ] = "Earth";

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, centralBodyMap );
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, simulationEndEpoch, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    Eigen::VectorXd finalIntegratedState = (--integrationResult.end( ) )->second;
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
