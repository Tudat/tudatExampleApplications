/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


/* A satellite is going to be launched on an elliptical orbit with min and max altitudes
respectively 180 km and 40000 km with inclination i=0 (and argument of perigee omega = 0).
A fixed target is set on the same plane at an altitude of 35000 km, fixed latitude of 30 deg.
What is the value of the RAAN for which the satellite achieves a minimum
approach distance from the target?*/

#include <pagmo/problem.hpp>
#include <pagmo/algorithms/sade.hpp>
#include <pagmo/algorithms/de1220.hpp>
#include <pagmo/algorithms/de.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>
#include "Problems/propagationTargeting.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

using namespace pagmo;

template< typename OutputStream, typename ScalarType,
          int NumberOfRows, int NumberOfColumns, int Options, int MaximumRows, int MaximumCols >
void writeValueToStream( OutputStream& stream, const Eigen::Matrix< ScalarType,
                         NumberOfRows, NumberOfColumns, Options,
                         MaximumRows, MaximumCols >& value,
                         const int precision, const std::string& delimiter,
                         const bool endLineAfterRow  = 0 )
{
    for ( int i = 0; i < value.rows( ); i++ )
    {
        for ( int j = 0; j < value.cols( ); j++ )
        {
            stream << delimiter << " "
                   << std::setprecision( precision ) << std::left
                   << std::setw( precision + 1 )
                   << value( i, j );


        }
        if( endLineAfterRow )
        {
            stream << std::endl;
        }
    }
    stream << std::endl;
}

template< typename ScalarType, int NumberOfRows, int NumberOfColumns >
void writeMatrixToFile( Eigen::Matrix< ScalarType, NumberOfRows, NumberOfColumns > matrixToWrite,
                        const std::string& outputFilename,
                        const int precisionOfMatrixEntries = 16,
                        const boost::filesystem::path& outputDirectory =
        "/home/dominic/Software/optimizationBundle/tudatBundle/tudatExampleApplications/libraryExamples/Pagmo2/bin/applications/",
                        const std::string& delimiter = "\t",
                        const std::string& header = "" )
{
    // Check if output directory exists; create it if it doesn't.
    if ( !boost::filesystem::exists( outputDirectory ) )
    {
        boost::filesystem::create_directories( outputDirectory );
    }

    // Open output file.
    std::string outputDirectoryAndFilename = outputDirectory.string( ) + "/" + outputFilename;
    std::ofstream outputFile_( outputDirectoryAndFilename.c_str( ) );

    // Write header
    outputFile_ << header;

    writeValueToStream( outputFile_, matrixToWrite, precisionOfMatrixEntries,
                        delimiter, true );

    outputFile_.close( );
}

void createGridSearch(
        pagmo::problem& problem,
        const std::vector< std::vector< double > >& bounds,
        const std::vector< int > numberOfPoints )
{
    Eigen::MatrixXd gridSearch = Eigen::MatrixXd( numberOfPoints.at( 0 ), 1 );

    double xSpacing = ( bounds[ 1 ][ 0 ] - bounds[ 0 ][ 0 ] ) / static_cast< double >( numberOfPoints.at( 0 ) - 1 );
    //double ySpacing = ( bounds[ 1 ][ 1 ] - bounds[ 0 ][ 1 ] ) / static_cast< double >( numberOfPoints.at( 1 ) - 1 );

    std::vector< double > decisionVector;
    decisionVector.resize( 1 );

    for( unsigned int i = 0; i < numberOfPoints.at( 0 ); i++ )
    {
        //for( unsigned int j = 0; j < numberOfPoints.at( 0 ); j++ )
        {
            decisionVector[ 0 ] = bounds[ 0 ][ 0 ] + static_cast< double >( i ) * xSpacing;
            //decisionVector[ 1 ] = bounds[ 0 ][ 1 ] + static_cast< double >( j ) * ySpacing;

            gridSearch( i, 0 ) = problem.fitness( decisionVector ).at( 0 );
            std::cout<<i<<" "<<gridSearch( i, 0 )<<std::endl;
        }
    }

    //writeMatrixToFile( gridSearch, "targetGridSearch.dat" );

}

int main()
{

    //Set seed for reproducible results
    pagmo::random_device::set_seed(255);

    //Load spice kernels
    std::string kernelsPath = tudat::input_output::getSpiceKernelPath( );
    tudat::spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    tudat::spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");

    // 1 - Instantiate a pagmo problem constructing it from a UDP
    // (user defined problem).
    problem prob{PropagationTargetingProblem( 180000.0, 40000000.0,
                                              35000000.0, 30.0 ) };

   // createGridSearch( prob, {{0.0}, {360.0}}, { 10000 } );

    // 2 - Instantiate a pagmo algorithm
    algorithm algo{de( )};

    // Create an island with 8 individuals
    island isl{algo, prob, 32};

    int i = 0;
    // For 30 generations optimise the population in the island
    for(  ; i < 500; i++ )
    {
         isl.evolve();
         while( isl.status()!=pagmo::evolve_status::idle )
             isl.wait();
         int c = isl.get_population().best_idx();
         vector_double cx = isl.get_population().champion_x();
         vector_double cf = isl.get_population().champion_f();
         print("GEN=", i, " ID=", c, " Separation=", cf[ 0 ], "m RAAN=", cx[ 0 ], "d\n" );
    }


//    // 3 - Instantiate an archipelago with 16 islands having each 20 individuals
//    archipelago archi{1, algo, prob, 16};

//    // 4 - Run the evolution in parallel on the 16 separate islands 10 times.
//    archi.evolve(4);
//    archi.wait();

//    // 6 - Print the decision and the fitness of the best solution in each island
//    for (const auto &isl : archi)
//    {
//        print(isl.get_population().champion_x(), "Selected RAAN (deg)\n");
//        print(isl.get_population().champion_f(), "Minimum separation from target (meters)\n");
//    }
}

