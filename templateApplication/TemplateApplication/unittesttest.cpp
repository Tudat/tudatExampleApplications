#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

int main( )
{


    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = physical_constants::JULIAN_YEAR;
    const double simulationEndEpoch = physical_constants::JULIAN_YEAR + 7.0 * tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    NamedBodyMap bodyMap = createBodies( bodySettings );

    double initialMass = 1200.0;

    // Create spacecraft object.
    bodyMap[ "LRO" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "LRO" ]->setConstantBodyMass( initialMass );

    bodyMap[ "LRO2" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "LRO2" ]->setConstantBodyMass( initialMass + 1000 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLRO;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLRO2;


    ///                                         ///
    ///  DEFINE YOUR ACCELERATION TYPES HERE    ///
    ///                                         ///
    ///                                         ///

    // Turn on central gravity of Sun
    accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                basic_astrodynamics::central_gravity ) );
    accelerationsOfLRO2[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
            basic_astrodynamics::central_gravity ) );

    // Turn on central gravity from Earth
    accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
            basic_astrodynamics::central_gravity ) );
    accelerationsOfLRO2[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
            basic_astrodynamics::central_gravity ) );

    // Turn on spherical harmonic perturbations of Moon
    accelerationsOfLRO[ "Moon" ].push_back(
            boost::make_shared< SphericalHarmonicAccelerationSettings >( 16, 16 ) );
    accelerationsOfLRO2[ "Moon" ].push_back(
            boost::make_shared< SphericalHarmonicAccelerationSettings >( 16, 16 ) );

    // Turn on central gravity of Mars
    accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
            basic_astrodynamics::central_gravity ) );
    accelerationsOfLRO2[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
            basic_astrodynamics::central_gravity ) );

    // Turn on central gravity of Jupiter
    accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
            basic_astrodynamics::central_gravity ) );
    accelerationsOfLRO2[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
            basic_astrodynamics::central_gravity ) );

    double thrustMagnitude = 1.0;
    double specificImpulse = 200.0;
    boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings =
               boost::make_shared< ConstantThrustEngineSettings >(
                        thrustMagnitude, specificImpulse );
    boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings2 =
               boost::make_shared< ConstantThrustEngineSettings >(
                        thrustMagnitude, specificImpulse );

    boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionSettings;
    boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionSettings2;

    thrustDirectionSettings = boost::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                "Moon", true, false );

    thrustDirectionSettings2 = boost::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                "Moon", true, false );

    accelerationsOfLRO[ "LRO" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
            thrustDirectionSettings, thrustMagnitudeSettings ) );
    accelerationsOfLRO2[ "LRO2" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
            thrustDirectionSettings2, thrustMagnitudeSettings2 ) );

    accelerationMap[ "LRO" ] = accelerationsOfLRO;
    accelerationMap[ "LRO2" ] = accelerationsOfLRO2;


    bodiesToPropagate.push_back( "LRO" );
    bodiesToPropagate.push_back( "LRO2" );

        ///                                                           ///
        ///  MODIFY CENTRAL BODY SETTINGS FOR SUB-QUESTION IF NEEDED  ///
        ///                                                           ///
        ///                                                           ///
    centralBodies.push_back( "Moon" );
    centralBodies.push_back( "Moon" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );



        // Set initial Kepler elements for LRO.
    Eigen::Vector6d initialKeplerElements;
    initialKeplerElements[ semiMajorAxisIndex ] = 1900.0E3;
    initialKeplerElements[ eccentricityIndex ] = 0.05;
    initialKeplerElements[ inclinationIndex ] =
            89.7 * mathematical_constants::PI / 180.0;


        ///                                                           ///
        ///  MODIFY INITIAL STATE ACCORDING TO STUDENT NUMBER         ///
        ///                                                           ///
        ///                                                           ///
        /// Student number: 4516893
        /// A=8, B=9, C=3
    initialKeplerElements[ argumentOfPeriapsisIndex ] = 36.0*8*mathematical_constants::PI/180.0;
    initialKeplerElements[ longitudeOfAscendingNodeIndex ] = 36.0*9*mathematical_constants::PI/180.0;
    initialKeplerElements[ trueAnomalyIndex ] = 36.0*3*mathematical_constants::PI/180.0;


    Eigen::Vector6d lroInitialCartesianState =  convertKeplerianToCartesianElements(
                initialKeplerElements, bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ) );


    initialKeplerElements[ trueAnomalyIndex ] = 0*3*mathematical_constants::PI/180.0;

    Eigen::Vector6d lroInitialCartesianState2 =  convertKeplerianToCartesianElements(
                initialKeplerElements, bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ) );


    Eigen::VectorXd initialCartesianState;
    initialCartesianState.resize(12);
    initialCartesianState.segment(0, 6) = lroInitialCartesianState;
    initialCartesianState.segment(6, 6) = lroInitialCartesianState2;


        ///                                                         ///
        ///  DEFINE PROPAGATION AND INTEGRATION SETTINGS HERE       ///
        ///                                                         ///
        ///                                                         ///
    double timeStep = 30.0;

    boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings;
    boost::shared_ptr< PropagatorSettings< double > > multiPropagatorSettings; // For question 1.6

    propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate,
                      initialCartesianState, simulationEndEpoch, cowell );

    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch );

    boost::shared_ptr< MassRateModelSettings > massRateModelSettings =
            boost::make_shared< FromThrustMassModelSettings >( 1 );
    boost::shared_ptr< MassRateModelSettings > massRateModelSettings2 =
            boost::make_shared< FromThrustMassModelSettings >( 1 );

    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;

    massRateModels[ "LRO" ] = createMassRateModel(
            "LRO", massRateModelSettings, bodyMap, accelerationModelMap );
    massRateModels[ "LRO2" ] = createMassRateModel(
            "LRO2", massRateModelSettings2, bodyMap, accelerationModelMap );

                // Create settings for propagating the mass of the vehicle.
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back( "LRO" );
    bodiesWithMassToPropagate.push_back( "LRO2" );

    Eigen::VectorXd initialBodyMasses = Eigen::VectorXd( 2 );
    initialBodyMasses( 0 ) = initialMass;
    initialBodyMasses( 1 ) = initialMass + 1000;

    boost::shared_ptr< PropagatorSettings< double > > massPropagatorSettings =
            boost::make_shared< MassPropagatorSettings< double > >(
                    bodiesWithMassToPropagate, massRateModels, initialBodyMasses, terminationSettings );

    std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;

    propagatorSettingsVector.push_back( massPropagatorSettings );
    propagatorSettingsVector.push_back( propagatorSettings );
    multiPropagatorSettings =
            boost::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsVector, terminationSettings );


    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, timeStep );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map< double, Eigen::VectorXd > cartesianIntegrationResult;

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
            bodyMap, integratorSettings, multiPropagatorSettings, true, false, false );

    cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );



    for( std::map< double, Eigen::VectorXd >::iterator it = cartesianIntegrationResult.begin();
         it != cartesianIntegrationResult.end(); ++it )
    {
        std::cout << it->second.transpose() << "\n";

    }

    return 0;

}
