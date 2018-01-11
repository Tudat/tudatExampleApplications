/*    Copyright (c) 2010-2017, Delft University of Technology
*    All rigths reserved
*
*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/

#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/moead.hpp"
#include "pagmo/algorithms/nsga2.hpp"
#include "pagmo/algorithms/ihs.hpp"
#include "Problems/earthMarsTransfer.h"


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

using namespace pagmo;

void createGridSearch(
        pagmo::problem& problem,
        const std::vector< std::vector< double > >& bounds,
        const std::vector< int > numberOfPoints )
{
    Eigen::MatrixXd gridSearch = Eigen::MatrixXd( numberOfPoints.at( 0 ), numberOfPoints.at( 1 ) );

    double xSpacing = ( bounds[ 1 ][ 0 ] - bounds[ 0 ][ 0 ] ) / static_cast< double >( numberOfPoints.at( 0 ) - 1 );
    double ySpacing = ( bounds[ 1 ][ 1 ] - bounds[ 0 ][ 1 ] ) / static_cast< double >( numberOfPoints.at( 1 ) - 1 );

    std::vector< double > decisionVector;
    decisionVector.resize( 2 );

    for( unsigned int i = 0; i < numberOfPoints.at( 0 ); i++ )
    {
        for( unsigned int j = 0; j < numberOfPoints.at( 0 ); j++ )
        {
            decisionVector[ 0 ] = bounds[ 0 ][ 0 ] + static_cast< double >( i ) * xSpacing;
            decisionVector[ 1 ] = bounds[ 0 ][ 1 ] + static_cast< double >( j ) * ySpacing;

            gridSearch( i, j ) = problem.fitness( decisionVector ).at( 0 );
        }
    }

    writeMatrixToFile( gridSearch, "porkchopPlot.dat" );

}

void printPopulationToFile( const int iterationIndex,
                            const std::vector< pagmo::vector_double >& population,
                            const bool isFitness )
{
    Eigen::MatrixXd matrixToPrint( population.size( ), population.at( 0 ).size( ) );
    for( unsigned int i = 0; i < population.size( ); i++ )
    {
        for( unsigned int j = 0; j < population.at( 0 ).size( ); j++ )
        {
            matrixToPrint( i, j ) = population.at( i ).at( j );
        }
    }

    if( !isFitness )
    {
        writeMatrixToFile( matrixToPrint, "population_lambert_mo_" + std::to_string( iterationIndex ) + ".dat" );
    }
    else
    {
        writeMatrixToFile( matrixToPrint, "fitness_lambert_mo_" + std::to_string( iterationIndex ) + ".dat" );
    }
}

//! Execute  main
int main( )
{
    // Set the PRNG seed, such that results are reproducable
    int seed = 123;
    pagmo::random_device::set_seed( seed );

    // We have two decision variables each with a lower and upper
    // bound, create a vector of vectors that will contain these.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 2, 0.0 ) );

    // Search between 2020 and 2025 for flight duration between 200
    // and 1000 days.
    bounds[ 0 ][ 0 ] = 2451545.0 - 10.0 * 365.0;
    bounds[ 1 ][ 0 ] = 2451545.0 + 10.0 * 365.0;
    bounds[ 0 ][ 1 ] = 100;
    bounds[ 1 ][ 1 ] = 2000;

    // Define the problem
    problem prob{EarthMarsTransfer( bounds, true )};

    //createGridSearch( prob, bounds, { 1000, 1000 } );
    // Select the self-adaptive differential evolution algorithm.
    // One generation per evolution step.
    algorithm algo{nsga2( )};

    // Create an island with 8 individuals
    //island isl{algo, prob, 16};
    island isl{algo, prob, 128};
    int i = 0;
    // For 30 generations optimise the population in the island
    for(  ; i < 100; i++ )
    {
        isl.evolve();

        printPopulationToFile( i, isl.get_population( ).get_x( ), false );
        printPopulationToFile( i, isl.get_population( ).get_f( ), true );

        //std::cout << "Best x: " << isl.get_population().champion_x()[0] << std::endl;
        //std::cout << "Best y: " << isl.get_population().champion_x()[1] << std::endl;
        while( isl.status()!=pagmo::evolve_status::idle )
            isl.wait();


    }

    return 0;

}
