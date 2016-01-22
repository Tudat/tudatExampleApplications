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

#include <string>
#include <iostream>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h>
#include <Tudat/External/SpiceInterface/spiceEphemeris.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

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
      std::map< double, tudat::basic_mathematics::Vector6d > stateHistory;
      double currentTime = 0.0;
      double endTime = (jdEnd - jdStart) * tudat::physical_constants::JULIAN_DAY;
      while( currentTime <= endTime )
      {
        // Get state and set in state history.
        stateHistory[ currentTime ] = ephem.getCartesianStateFromEphemeris(
            currentTime, jdStart );

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
