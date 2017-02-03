/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include <string>
#include <iostream>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h>
#include <Tudat/External/SpiceInterface/spiceEphemeris.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Basics/basicTypedefs.h>

#include <json/json.h>
#include <json/value.h>

//! Execute propagation of orbits of Asterix and Obelix around the Earth.
int main( )
{
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::ephemerides;
    using namespace tudat::input_output;
    using namespace tudat::spice_interface;

    // Test JSON
    // __FILE__ only gives the absolute path in the header file!
    std::string appPath( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    appPath = appPath.substr( 0, appPath.find_last_of("/\\")+1);

    std::cout << "Reading config file:" << std::endl
              << "\t" << appPath + "config.json" << std::endl;

    // Construct Json Database
    Json::Value jsonRoot;

    std::ifstream config_doc(appPath + "config.json", std::ifstream::binary);
    config_doc >> jsonRoot; // write file to jsonRoot


    Json::Value bodies = jsonRoot["bodies"];
    Json::Value gdStart = jsonRoot["sim"]["start"];
    double jdStart = convertCalendarDateToJulianDay(gdStart[0].asInt(), gdStart[1].asInt(),
                                                    gdStart[2].asInt(), gdStart[3].asInt(),
                                                    gdStart[4].asInt(), gdStart[5].asInt());
    double timeStep = jsonRoot["sim"]["step"].asDouble();
    Json::Value gdEnd = jsonRoot["sim"]["end"];
    double jdEnd = convertCalendarDateToJulianDay(gdEnd[0].asInt(), gdEnd[1].asInt(), gdEnd[2].asInt(),
                                                  gdEnd[3].asInt(), gdEnd[4].asInt(), gdEnd[5].asInt());
    std::string refFrame(jsonRoot["defs"].get("frame", "ICRF").asString());
    std::string refEpoch(jsonRoot["defs"].get("epoch", "J2000").asString());
    std::string refObser(jsonRoot["defs"].get("observer", "SSB").asString());

    // Do some printing
    std::cout << "Bodies:   ";
    for ( unsigned int i = 0; i < bodies.size(); i++ ){
      std::cout << bodies[i].asString() << ", ";
    }
    std:: cout << std::endl;
    //std::cout << "Bodies: " << jsonRoot["bodies"].asString() << std::endl;
    std::cout << "Frame:    " << refFrame << std::endl;
    std::cout << "Observer: " << refObser << std::endl;
    std::cout << "Epoch:    " << refEpoch << std::endl;
    std::cout << "Start:    " << jdStart << std::endl;
    std::cout << "End:      " << jdEnd << std::endl;
    std::cout << "Diff:     " << (jdEnd - jdStart) << std::endl;
    std::cout << "Step:     " << timeStep << std::endl;

    // Load Spice kernels.
    loadSpiceKernelInTudat( getTudatRootPath( ) + "External/SpiceInterface/Kernels/de421.bsp" );
    loadSpiceKernelInTudat( getTudatRootPath( ) + "External/SpiceInterface/Kernels/pck00009.tpc" );
    loadSpiceKernelInTudat( getTudatRootPath( ) + "External/SpiceInterface/Kernels/naif0009.tls" );
    loadSpiceKernelInTudat( getTudatRootPath( ) +
                            "External/SpiceInterface/Kernels/de-403-masses.tpc" );

    for ( unsigned int i = 0; i < bodies.size(); i++ ){
      std::string bodyName = bodies[i].asString();
      SpiceEphemeris ephem(bodyName, refObser, false, false, false, refEpoch );

      std::cout << "Compute ephemeri" << std::endl;
      std::map< double, Eigen::Vector6d > stateHistory;
      double currentTime = 0.0;
      double endTime = (jdEnd - jdStart) * tudat::physical_constants::JULIAN_DAY;
      while( currentTime <= endTime )
      {
        // Get state and set in state history.
        stateHistory[ currentTime ] = ephem.getCartesianState(
            currentTime + ( jdStart - tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                    tudat::physical_constants::JULIAN_DAY );

        // Update current time.
        currentTime += timeStep;
      }


      std::cout << "Writing out" << std::endl
                << "\t" << appPath + bodyName + "History.dat" << std::endl;
      writeDataMapToTextFile( stateHistory, bodyName + "History.dat", appPath,
                              "", 10, 10, "\t");
    }
    return 0;
}
