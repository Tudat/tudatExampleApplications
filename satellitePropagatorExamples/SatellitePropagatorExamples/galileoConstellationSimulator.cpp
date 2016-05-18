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
 *      130250    K. Kumar          Updated gravitational acceleration model references; renamed
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
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

#include <Tudat/Astrodynamics/Propagators/dynamicsSimulator.h>
#include <Tudat/External/SpiceInterface/spiceInterface.h>
#include <Tudat/SimulationSetup/body.h>
#include <Tudat/SimulationSetup/createAccelerationModels.h>
#include <Tudat/SimulationSetup/defaultBodies.h>

#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

//! Execute simulation of Galileo constellation around the Earth.
int main( )
{
    using namespace tudat;

    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;

    using orbital_element_conversions::xCartesianPositionIndex;
    using orbital_element_conversions::yCartesianPositionIndex;
    using orbital_element_conversions::zCartesianPositionIndex;
    using orbital_element_conversions::xCartesianVelocityIndex;
    using orbital_element_conversions::yCartesianVelocityIndex;
    using orbital_element_conversions::zCartesianVelocityIndex;

    using basic_mathematics::Vector6d;

    using gravitation::CentralJ2J3J4GravitationalAccelerationModel;

    using input_output::DoubleKeyTypeVectorXdValueTypeMap;
    using input_output::writeDataMapToTextFile;

    using orbital_element_conversions::convertKeplerianToCartesianElements;

    using numerical_integrators::RungeKutta4Integrator;

    using unit_conversions::convertDegreesToRadians;



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
    const double simulationEndEpoch = 1.0 * physical_constants::JULIAN_DAY;

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
    const double earthJ3 = -0.0000050323;
    const double earthJ4 = -0.0000016204;
    const double earthEquatorialRadius = 6378.1363e3;

    // Set orbital parameters of Galileo constellation.
    const double semiMajorAxis = 23222.0e3 + 6378.1e3;                                // [km]
    const double eccentricity = 0.0;                                                  // [-]
    const double inclination = convertDegreesToRadians( 56.0 );                       // [rad]
    const double argumentOfPeriapsis = 0.0;                                           // [rad]
    const double longitudeOfAscendingNodeSpacing
            = 2.0 * mathematical_constants::PI / numberOfPlanes;                          // [rad]
    const double trueAnomalySpacing
            = 2.0 * mathematical_constants::PI / numberOfSatellitesPerPlane;              // [rad]

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

    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Create list of satellites and acceleration models per satellite.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize, simulationEndEpoch + 10.0 * fixedStepSize );
    Eigen::Matrix< double, 50, 50 > earthCosineCoefficients = Eigen::Matrix< double, 50, 50 >::Zero( );
    Eigen::Matrix< double, 50, 50 > earthSineCoefficients = Eigen::Matrix< double, 50, 50 >::Zero( );
    earthCosineCoefficients( 0, 0 ) = 1.0;
    earthCosineCoefficients << 1.0, 0.0,
            ( 1.0 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) * earthJ2 ),
            ( 1.0 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 0 ) * earthJ3 ),
            ( 1.0 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 4, 0 ) * earthJ4 );

    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
                basic_mathematics::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< simulation_setup::SphericalHarmonicsGravityFieldSettings >(
                earthGravitationalParameter, earthEquatorialRadius, earthCosineCoefficients, earthSineCoefficients, "IAU_Earth" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = NULL;
    bodySettings[ "Earth" ]->shapeModelSettings = NULL;

    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    SelectedAccelerationMap accelerationMap;

    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralBodyMap;
    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6 * numberOfSatellites );

    std::string currentSatelliteName;
    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        currentSatelliteName =  "Satellite" + boost::lexical_cast< std::string >( i );
        bodyMap[ currentSatelliteName ] = boost::make_shared< simulation_setup::Body >( );

        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfCurrentSatellite;
        accelerationsOfCurrentSatellite[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
        accelerationMap[ currentSatelliteName ] = accelerationsOfCurrentSatellite;

        bodiesToPropagate.push_back( currentSatelliteName );
        centralBodies.push_back( "Earth" );
        centralBodyMap[ currentSatelliteName ] = "Earth";

        systemInitialState.segment( i * 6, 6 ) = initialConditions.col( i );
    }

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
//    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

//    std::vector< std::map< double, Eigen::VectorXd > > allSatellitesPropagationHistory;
//    allSatellitesPropagationHistory.resize( numberOfSatellites );

//    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
//         stateIterator != integrationResult.end( ); stateIterator++ )
//    {
//        for( unsigned int i = 0; i < allSatellitesPropagationHistory.size( ); i++ )
//        {
//            allSatellitesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
//        }
//    }

//    // Write results to file.

//    // Loop over all satellites.
//    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
//    {
//        // Set filename for output data.
//        std::stringstream outputFilename;
//        outputFilename << "galileoSatellite" << i + 1 << ".dat";

//        // Write propagation history to file.
//        writeDataMapToTextFile( allSatellitesPropagationHistory.at( i ),
//                                outputFilename.str( ),
//                                outputDirectory,
//                                "",
//                                std::numeric_limits< double >::digits10,
//                                std::numeric_limits< double >::digits10,
//                                ", " );
//    }

    ///////////////////////////////////////////////////////////////////////////

    return 0;
}
