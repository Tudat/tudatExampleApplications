/*    Copyright (c) 2010-2018, Delft University of Technology
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
 *      120522    A. Ronse          First creation of code.
 *
 *    References
 *      Williams, Dr. David R., "Moon Fact Sheet", NASA (National Space Science Data Center),
 *         http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html, last accessed: 22 May 2012
 *
 *    Notes
 *
 */

#include <ctime>
#include <iostream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>

int main( )
{
    // Use Boost library to make a random number generator.
    boost::mt19937 randomNumbergenerator( time( 0 ) );
    boost::random::uniform_real_distribution< > uniformDistribution( 0.0, 10.0 );
    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
            generateRandomNumbers( randomNumbergenerator, uniformDistribution );

    // Generate random altitude value between 0 and 10 km.
    const double altitudeKilometers = generateRandomNumbers( );

    // Use the Tudat Core library to convert km to m.
    const double altitudeMeters = tudat::unit_conversions::convertKilometersToMeters(
            altitudeKilometers );

    // Use the Eigen library to create position vectors.
    Eigen::Vector3d positionOfBodySubjectToAcceleration;
    positionOfBodySubjectToAcceleration << 1737.1e3 + altitudeMeters, 0.0, 0.0;
    const Eigen::Vector3d positionOfBodyExertingAcceleration = Eigen::Vector3d::Zero( );

    // Use the Tudat library to compute the acceleration vector.
    const Eigen::Vector3d gravitationalAcceleration =
            tudat::gravitation::computeGravitationalAcceleration(
                    positionOfBodySubjectToAcceleration, 4.9e12,
                    positionOfBodyExertingAcceleration );

    // Print the altitude and norm of the acceleration vector.
    std::cout << "Hello world!" << std::endl;
    std::cout << "I am floating " << altitudeKilometers << " km above the Moon's surface." <<
            std::endl;
    std::cout << "The gravitational acceleration here is " <<
            gravitationalAcceleration.norm() << " m/s^2."  << std::endl;

    return 0;
}
