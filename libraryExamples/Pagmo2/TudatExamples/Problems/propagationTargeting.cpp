/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include"propagationTargeting.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>


PropagationTargetingProblem::PropagationTargetingProblem( const double altitudeOfPerigee,
        const double altitudeOfApogee, const double altitudeOfTarget, const double longitudeOfTarget ) :
    altitudeOfPerigee_( altitudeOfPerigee ), altitudeOfApogee_( altitudeOfApogee ),
    altitudeOfTarget_( altitudeOfTarget ), longitudeOfTarget_( longitudeOfTarget ) { }


std::vector<double> PropagationTargetingProblem::fitness(const std::vector<double> &x) const
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::input_output;

    // Definition of the orbit
    const double earthRadius = spice_interface::getAverageRadius( "Earth" );
    const double radiusOfPerigee =  earthRadius + altitudeOfPerigee_;
    const double radiusOfApogee = earthRadius + altitudeOfApogee_;
    const double earthGravitationalParameter = spice_interface::getBodyGravitationalParameter( "Earth" );
    const double semiMajorAxis = (radiusOfApogee + radiusOfPerigee)/2.0;

    //Integration time: half a orbit
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = mathematical_constants::PI *
            sqrt(pow(semiMajorAxis,3)/earthGravitationalParameter);
    const double fixedStepSize = 30;

    // Create the body Earth from Spice interface
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( {"Earth"} );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = NULL;
    bodySettings[ "Earth" ]->shapeModelSettings = NULL;
    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "J2000" );

    //Create bodyMap and add the satellite as an empty body
    NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );
    bodyMap["Satellite"] = boost::make_shared<Body>();

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    //Define position of the target at 35000 km from Earth at 30 deg latitude
    Eigen::Vector3d target;
    target[0] = (earthRadius + altitudeOfTarget_) * cos(longitudeOfTarget_*mathematical_constants::PI/180);
    target[1] = (earthRadius + altitudeOfTarget_) * sin(longitudeOfTarget_*mathematical_constants::PI/180);
    target[2] = 0.0;

    //Define initial position of satellite at the perigee
    Eigen::Vector6d initialKeplerElements;
    initialKeplerElements[ semiMajorAxisIndex ] = semiMajorAxis;
    initialKeplerElements[ eccentricityIndex ] = (radiusOfApogee - radiusOfPerigee)/(radiusOfApogee + radiusOfPerigee);
    initialKeplerElements[ inclinationIndex ] = 0.0;
    initialKeplerElements[ argumentOfPeriapsisIndex ] = 0.0;
    initialKeplerElements[ longitudeOfAscendingNodeIndex ] = x[0] * mathematical_constants::PI/180.0;
    initialKeplerElements[ trueAnomalyIndex ] = 0.0; //start at perigee

    const Eigen::Vector6d systemInitialState = convertKeplerianToCartesianElements(
                initialKeplerElements, earthGravitationalParameter );

    //Setup simulation. Simple Keplerian orbit (only central-gravity of Earth)
    std::vector< std::string > bodiesToPropagate = { "Satellite" };
    std::vector< std::string > centralBodies = { "Earth" };
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    accelerationsOfSatellite[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                              basic_astrodynamics::central_gravity ) );
    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    //Setup propagator (cowell) and integrator (RK4 fixed stepsize)
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    //Start simulation
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );

    //Retrieve results
    std::map< double, Eigen::VectorXd > integrationResult =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    //Find minimum distance from target
    Eigen::Vector3d separationFromTarget = integrationResult.begin( )->second.segment( 0, 3 ) - target;
    double bestDistanceFromTarget = sqrt(pow(separationFromTarget[0],2) + pow(separationFromTarget[1],2) +
            + pow(separationFromTarget[3],2));
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
         stateIterator != integrationResult.end( ); stateIterator++ )
    {
        separationFromTarget = stateIterator->second.segment( 0, 3 ) - target;
        const double distanceFromTarget = separationFromTarget.norm( );

        if( distanceFromTarget < bestDistanceFromTarget )
        {
            bestDistanceFromTarget = distanceFromTarget;
        }

    }
    std::vector< double > output = {bestDistanceFromTarget} ;
    return output;


}


std::pair<std::vector<double>, std::vector<double>> PropagationTargetingProblem::get_bounds() const
{
    return {{0.0}, {360.0}};
}
