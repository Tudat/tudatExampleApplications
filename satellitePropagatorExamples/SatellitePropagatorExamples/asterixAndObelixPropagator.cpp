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
 *      130107    K. Kumar          Updated license in file header.
 *      130225    K. Kumar          Updated gravitational acceleration model references; renamed
 *                                  file; fixed error in assigning Obelix state derivative model;
 *                                  made variables const-correct.
 *
 *    References
 *
 *    Notes
 *
 */

#include <boost/assign/list_of.hpp>
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
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include <Tudat/Astrodynamics/Propagators/dynamicsSimulator.h>
#include <Tudat/External/SpiceInterface/spiceInterface.h>
#include <Tudat/SimulationSetup/body.h>
#include <Tudat/SimulationSetup/createAccelerationModels.h>
#include <Tudat/SimulationSetup/defaultBodies.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <utility>

#include <Eigen/Core>

//! Execute propagation of orbits of Asterix and Obelix around the Earth.
int main( )
{
    using namespace tudat;

    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;

    using tudat::basic_astrodynamics::AccelerationModel3dPointer;
    using tudat::orbital_element_conversions::semiMajorAxisIndex;
    using tudat::orbital_element_conversions::eccentricityIndex;
    using tudat::orbital_element_conversions::inclinationIndex;
    using tudat::orbital_element_conversions::argumentOfPeriapsisIndex;
    using tudat::orbital_element_conversions::longitudeOfAscendingNodeIndex;
    using tudat::orbital_element_conversions::trueAnomalyIndex;
    using tudat::orbital_element_conversions::xCartesianPositionIndex;
    using tudat::orbital_element_conversions::yCartesianPositionIndex;
    using tudat::orbital_element_conversions::zCartesianPositionIndex;
    using tudat::orbital_element_conversions::xCartesianVelocityIndex;
    using tudat::orbital_element_conversions::yCartesianVelocityIndex;
    using tudat::orbital_element_conversions::zCartesianVelocityIndex;

    using tudat::basic_mathematics::Vector6d;

    using tudat::input_output::writeDataMapToTextFile;

    using tudat::numerical_integrators::RungeKutta4Integrator;

    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;

    using tudat::unit_conversions::convertDegreesToRadians;

    ///////////////////////////////////////////////////////////////////////////

    // Input deck.

    // Set output directory.
    const std::string outputDirectory = "/home/dominicdirkx/Software/JaccoNewBundle/tudatBundle/tudatExampleApplications/satellitePropagatorExamples/bin/applications";
    if( outputDirectory == "" )
    {
        std::cerr<<"Error, output directory not specified (modify outputDirectory variable to "<<
                   " required output directory)."<<std::endl;
    }

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 60.0;

    // Set initial conditions for satellites that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to
    // Cartesian elements.

    // Set Keplerian elements for Asterix.
    Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0e3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Set Keplerian elements for Obelix.
    Vector6d obelixInitialStateInKeplerianElements( 6 );
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
    const Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );

    // Convert Obelix state from Keplerian elements to Cartesian elements.
    const Vector6d obelixInitialState = convertKeplerianToCartesianElements(
                obelixInitialStateInKeplerianElements,
                earthGravitationalParameter );

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////


    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Create list of satellites and acceleration models per satellite.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize, simulationEndEpoch + 10.0 * fixedStepSize );
    Eigen::Matrix< double, 5, 1 > earthCosineCoefficients;
    earthCosineCoefficients << 1.0, 0.0,
            ( 1.0 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) * earthJ2 ),
            ( 1.0 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 0 ) * earthJ3 ),
            ( 1.0 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 4, 0 ) * earthJ4 );

    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
                basic_mathematics::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< simulation_setup::SphericalHarmonicsGravityFieldSettings >(
                earthGravitationalParameter, earthEquatorialRadius, earthCosineCoefficients, Eigen::Matrix< double, 5, 1 >::Zero( ), "IAU_Earth" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = NULL;
    bodySettings[ "Earth" ]->shapeModelSettings = NULL;

    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    SelectedAccelerationMap accelerationMap;

    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralBodyMap;
    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 12 );

    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Obelix" ] = boost::make_shared< simulation_setup::Body >( );

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationMap[  "Asterix" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );
    centralBodyMap[  "Asterix" ] = "Earth";

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfObelix;
    accelerationsOfObelix[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationMap[  "Obelix" ] = accelerationsOfObelix;
    bodiesToPropagate.push_back( "Obelix" );
    centralBodies.push_back( "Earth" );
    centralBodyMap[  "Obelix" ] = "Earth";

    systemInitialState.segment( 0, 6 ) = asterixInitialState;
    systemInitialState.segment( 6, 6 ) = obelixInitialState;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, centralBodyMap );
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, simulationEndEpoch, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    std::map< double, Eigen::VectorXd > asterixPropagationHistory;
    std::map< double, Eigen::VectorXd > obelixPropagationHistory;

    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
         stateIterator != integrationResult.end( ); stateIterator++ )
    {
        asterixPropagationHistory[ stateIterator->first ] = stateIterator->second.segment( 0, 6 );
        obelixPropagationHistory[ stateIterator->first ] = stateIterator->second.segment( 6, 6 );
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
