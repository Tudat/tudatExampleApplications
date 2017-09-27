/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/JsonInterface/jsonInterface.h>

#include <Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h>
#include <SatellitePropagatorExamples/applicationOutput.h>

class ApolloJsonSimulationManager : public tudat::json_interface::JsonSimulationManager< >
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
        getBody( "Apollo" )->setAerodynamicCoefficientInterface( tudat::unit_tests::getApolloCoefficientInterface( ) );
    }

    // Override resetExportSettings method
    virtual void resetExportSettings( )
    {
        // First, call the original resetExportSettings, which uses the information in the JSON file
        JsonSimulationManager::resetExportSettings( );

        // Then, replace the output file paths (empty strings placeholders had been specified in the JSON file)
        const std::string outputDirectory = tudat_applications::getOutputPath( ) + "ApolloCapsuleExampleJSON/";
        getExportSettings( 0 )->setOutputFile( outputDirectory + "apolloPropagationHistory.dat" );
        getExportSettings( 1 )->setOutputFile( outputDirectory + "apolloDependentVariableHistory.dat" );
    }

    // Override resetPropagatorSettings method
    virtual void resetPropagatorSettings( )
    {
        // First, call the original resetPropagatorSettings, which uses the information in the JSON file
        JsonSimulationManager::resetPropagatorSettings( );

        // Define constant 30 degree angle of attack
        double constantAngleOfAttack = 30.0 * tudat::mathematical_constants::PI / 180.0;
        getBody( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->
                setOrientationAngleFunctions( boost::lambda::constant( constantAngleOfAttack ) );

        // Set spherical elements for Apollo.
        using namespace tudat::orbital_element_conversions;
        Eigen::Vector6d apolloSphericalEntryState;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
                getBody( "Earth" )->getShapeModel( )->getAverageRadius( ) + 120.0E3;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.0;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.7E3;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
                -0.9 * tudat::mathematical_constants::PI / 180.0;
        apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

        // Convert apollo state from spherical elements to Cartesian elements.
        Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState( apolloSphericalEntryState );
        systemInitialState = tudat::ephemerides::transformStateToGlobalFrame(
                    systemInitialState, getStartEpoch( ), getBody( "Earth" )->getRotationalEphemeris( ) );

        // Reset initial states (zero placeholder vector had been specified in the JSON file)
        getPropagatorSettings( )->resetInitialStates( systemInitialState );
    }
};

//! Execute propagation of orbits of Apollo during entry using the JSON Interface.
int main( )
{
    const std::string cppFilePath( __FILE__ );
    const std::string cppFolder = cppFilePath.substr( 0, cppFilePath.find_last_of("/\\") + 1 );

    using namespace tudat::json_interface;

    ApolloJsonSimulationManager jsonSimulationManager( cppFolder + "apolloCapsuleEntry.json" );
    jsonSimulationManager.updateSettingsFromJsonObject( );
    jsonSimulationManager.run( );
    jsonSimulationManager.exportResults( );

    return EXIT_SUCCESS;
}

