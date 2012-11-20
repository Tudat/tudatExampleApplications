/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120522    A. Ronse          First creation of code.
 *
 *    References
 *      Williams, Dr. David R., "Moon Fact Sheet", NASA (National Space Science Data Center),
 *         http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html, last accessed: 22 May 2012
 *
 */

#include <iostream>
#include <ctime>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

int main( )
{
    // Use Boost library to make a random number generator.
    boost::mt19937 randomNumbergenerator( time( 0 ) );
    boost::random::uniform_real_distribution< > uniformDistribution( 0.0, 10.0 );
    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
            generateRandomNumbers( randomNumbergenerator, uniformDistribution );

    // Generate random altitude value between 0 and 10 km.
    double altitudeKilometers = generateRandomNumbers( );

    // Use the TudatCore library to convert Km to m
    double altitudeMeters = tudat::unit_conversions::convertKilometersToMeters(
            altitudeKilometers );

    // Use the Eigen library to create position vectors.
    Eigen::Vector3d positionOfBodySubjectToAcceleration;
    positionOfBodySubjectToAcceleration << 1737.1e3 + altitudeMeters, 0, 0;
    Eigen::Vector3d positionOfBodyExertingAcceleration = Eigen::Vector3d::Zero( );

    // Use the Tudat library to compute the acceleration vector.
    Eigen::Vector3d gravitationalAcceleration =
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
