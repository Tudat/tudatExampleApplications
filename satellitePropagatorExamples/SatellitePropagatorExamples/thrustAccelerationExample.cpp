// This is the main file
//#include <Eigen/Core>
//#include <iostream>
//#include <cmath>

//#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelElementConversions.h"
//#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
//#include "Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h"
//#include <Tudat/Basics/testMacros.h>
//#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
//#include <Tudat/External/SpiceInterface/spiceEphemeris.h>
//#include <Tudat/External/SpiceInterface/spiceRotationalEphemeris.h>
//#include <Tudat/InputOutput/basicInputOutput.h>
//#include <Tudat/SimulationSetup/EnvironmentSetup/body.h>
//#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
//#include <Tudat/SimulationSetup/PropagationSetup/createMassRateModels.h>
//#include <Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h>
//#include <Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h>
//#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/applicationOutput.h>

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>


int main()
{


    // FOURTH TRY

    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace unit_conversions;

    /// CREATE ENVIRONMENT AND VEHICLE

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 60.0;

    // Define body settings for simulation.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                basic_mathematics::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );


    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    /// CREATE ACCELERATIONS

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Find filepath of this cpp file
    std::string cppFilePath( __FILE__ );

    // Find folder of this cpp file
    std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    Eigen::MatrixXd thrustForceMatrix =
            input_output::readMatrixFromFile( cppFolder + "unitTests/nonconstantThrust.txt" , " \t", "#" );

    // Make map for thrust data
    std::map< double, Eigen::Matrix3d > thrustData;

    // Fill thrustData map using thrustForceMatrix Eigen matrix
    for ( int i = 0; i < thrustForceMatrix.rows(); i++ )
    {
        Eigen::Matrix3d blabla = thrustForceMatrix.rightCols( 3 ).row( i );
        thrustData[ thrustForceMatrix( i, 0 ) ] = blabla;
    }


    // Define order of Lagrange interpolator
    int interpolatorOrder = 3;

    // Make interpolator
    boost::shared_ptr< interpolators::LagrangeInterpolatorSettings >
            thrustInterpolatorSettingsPointer = boost::make_shared< interpolators::LagrangeInterpolatorSettings >
            ( interpolatorOrder, false, interpolators::huntingAlgorithm,
              interpolators::lagrange_cubic_spline_boundary_interpolation );

    //!
    // Creating settings for thrust force
    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > >
            thrustInterpolatorPointer = interpolators::createOneDimensionalInterpolator< double, Eigen::Vector3d >(
                thrustData, thrustInterpolatorSettingsPointer );

    boost::function< double( const double ) > specificImpulseFunction = boost::lambda::constant( 3000.0 ); //? See function


    //!
    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Vehicle" ].push_back(
                boost::make_shared< tudat::simulation_setup::ThrustAccelerationSettings >(
                    thrustInterpolatorPointer,
                    specificImpulseFunction ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    /// CREATE PROPAGATION SETTINGS

    // Set initial conditions for the vehicle satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.

    // Set Keplerian elements for vehicle.
    Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 72130.0e3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.6;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 169.0 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 45.0 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 80.0 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 15.0 );

    double vehicleMass = 2000.0;

    // Convert vehicle state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                vehicleInitialStateInKeplerianElements,
                earthGravitationalParameter );

    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = (
                createMassRateModel( "Vehicle", boost::make_shared< tudat::simulation_setup::FromThrustMassModelSettings >( 1 ),
                                     bodyMap, accelerationModelMap ) );
    // Is the above expression specifying how the mass should be propagated?

    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                    boost::make_shared< PropagationTimeTerminationSettings >( 1000.0 );
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, fixedStepSize );

    boost::shared_ptr< PropagatorSettings< double > > massPropagatorSettings =
            boost::make_shared< MassPropagatorSettings< double > >(
                boost::assign::list_of( "Vehicle" ), massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

    std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

    /// PROPAGATE ORBIT

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );


    /// PROVIDE OUTPUT TO CONSOLE AND FILES

    Eigen::VectorXd finalIntegratedState = (--integrationResult.end( ) )->second;
    // Print the position (in km) and the velocity (in km/s) at t = 0.
    std::cout << "Single Earth-Orbiting Satellite with subject to varying thrust force Example." << std::endl <<
                 "The initial position vector of Vehicle is [km]:" << std::endl <<
                 systemInitialState.segment( 0, 3 ) / 1E3 << std::endl <<
                 "The initial velocity vector of Vehicle is [km/s]:" << std::endl <<
                 systemInitialState.segment( 3, 3 ) / 1E3 << std::endl;

    // Print the position (in km) and the velocity (in km/s) at t = 86400.
    std::cout << "After " << simulationEndEpoch <<
                 " seconds, the position vector of Vehicle is [km]:" << std::endl <<
                 finalIntegratedState.segment( 0, 3 ) / 1E3 << std::endl <<
                 "And the velocity vector of Vehicle is [km/s]:" << std::endl <<
                 finalIntegratedState.segment( 3, 3 ) / 1E3 << std::endl;
}
