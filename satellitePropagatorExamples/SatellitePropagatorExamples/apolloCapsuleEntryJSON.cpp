/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/InputOutput/JsonInterface/simulation.h>

#include <Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h>

#include <SatellitePropagatorExamples/applicationOutput.h>

namespace tudat
{

namespace json_interface
{

class ApolloJsonSimulationManager : public JsonSimulationManager< >
{
public:
    // Inherit constructor.
    using JsonSimulationManager< >::JsonSimulationManager;

protected:
    // Override resetBodies method
    virtual void resetBodies( )
    {
        // First, call the original resetBodies, which uses the information in the JSON file
        JsonSimulationManager::resetBodies( );

        // Then, create vehicle's aerodynamic coefficients interface
        getBody( "Apollo" )->setAerodynamicCoefficientInterface( unit_tests::getApolloCoefficientInterface( ) );
    }

    // Override resetPropagatorSettings method
    virtual void resetPropagatorSettings( )
    {
        // First, call the original resetPropagatorSettings, which uses the information in the JSON file
        JsonSimulationManager::resetPropagatorSettings( );

        using namespace tudat::orbital_element_conversions;

        // Define constant 30 degree angle of attack
        double constantAngleOfAttack = 30.0 * mathematical_constants::PI / 180.0;
        getBody( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->
                setOrientationAngleFunctions( boost::lambda::constant( constantAngleOfAttack ) );

        // Set spherical elements for Apollo.
        Eigen::Vector6d apolloSphericalEntryState;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
                spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.0;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.7E3;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
                -0.9 * mathematical_constants::PI / 180.0;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

        // Convert apollo state from spherical elements to Cartesian elements.
        Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState( apolloSphericalEntryState );
        systemInitialState = ephemerides::transformStateToGlobalFrame(
                    systemInitialState, getStartEpoch( ), getBody( "Earth" )->getRotationalEphemeris( ) );

        // Reset initial states (zero placeholder vector had been specified in the JSON file)
        getPropagatorSettings( )->resetInitialStates( systemInitialState );
    }
};

}  // namespace json_interface

}  // namespace tudat


//! Execute propagation of orbits of Apollo during entry using the JSON Interface.
int main( )
{
    using namespace tudat::json_interface;

    const std::string cppFilePath( __FILE__ );
    const std::string cppFolder = cppFilePath.substr( 0, cppFilePath.find_last_of("/\\") + 1 );

    ApolloJsonSimulationManager jsonSimulationManager;
    jsonSimulationManager.setUpFromJSONFile( cppFolder + "apolloCapsuleEntry.json" );
    jsonSimulationManager.run( );
    jsonSimulationManager.exportResults( );

    return EXIT_SUCCESS;
}

