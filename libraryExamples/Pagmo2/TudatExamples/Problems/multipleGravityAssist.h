/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EXAMPLE_PAGMO_MULTIPLE_GRAVITY_ASSIST_H
#define TUDAT_EXAMPLE_PAGMO_MULTIPLE_GRAVITY_ASSIST_H

//#include <vector>
//#include <utility>
//#include <limits>

//#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
//#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
//#include <Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h>

//#include "mga.h"

//#include "pagmo/island.hpp"
//#include "pagmo/io.hpp"
//#include "pagmo/serialization.hpp"
//#include "pagmo/problem.hpp"

//#include <Eigen/Core>

//using namespace pagmo;

////! Test function for a new interplanetary trajectory class in Tudat
//struct MultipleGravityAssist
//{
//    MultipleGravityAssist(
//            const std::vector< int >& flybySequence,
//            const std::pair< std::vector< double >, std::vector< double > >& problemBounds )
//    {
//        mgaObject_.type = total_DV_rndv;
//        mgaObject_.sequence = flybySequence;
//        mgaObject_.rev_flag.resize( flybySequence.size( ) );

//        for( unsigned int i = 0; i < flybySequence.size( ); i++ )
//        {
//            mgaObject_.rev_flag[ i ] = 0;
//        }

//        mgaObject_.Isp = 300.0;
//        mgaObject_.mass = 8000.0;
//        mgaObject_.DVlaunch = 0.0;

////                double e;							//insertion e (only in case total_DV_orbit_insertion)
////                double rp;							//insertion rp in km (only in case total_DV_orbit_insertion)
////                customobject asteroid;				//asteroid data (in case fly-by sequence has a final number = 10)
//    }

//    // Calculates the fitness
//    std::vector< double > fitness( const std::vector< double > &x ) const
//    {
//        std::vector<double> rp;
//        std::vector<double> DV;
//        double obj_funct;
//        int result = MGA( x, mgaObject_, rp, DV, obj_funct );
//        return { obj_funct };
//    }

//    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
//    {
//        return problemBounds_;
//    }

//    std::string get_name( ) const
//    {
//        return "Multiple Gravity Assist";
//    }

//    template <typename Archive>
//    void serialize(Archive &ar)
//    {
//        ar( problemBounds_);
//    }

//    vector_double::size_type get_nobj() const
//    {
//        return 1u;
//    }

//private:

//    mgaproblem mgaObject_;

//    const std::pair< std::vector< double >, std::vector< double > > problemBounds_;
//};

#include <vector>
#include <utility>
#include <limits>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include <Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h>

#include "mga.h"

#include "pagmo/island.hpp"
//#include "pagmo/external/cereal/cereal.hpp"
#include "pagmo/io.hpp"
#include "pagmo/serialization.hpp"
#include "pagmo/problem.hpp"

#include <Eigen/Core>

typedef Eigen::Matrix< double, 6, 1 > StateType;

using namespace pagmo;

//! Test function for a new interplanetary trajectory class in Tudat
struct MultipleGravityAssist
{

    MultipleGravityAssist( const bool useTripTime = false ): useTripTime_( useTripTime ){ }

    MultipleGravityAssist( std::vector< std::vector< double > > &bounds,
                           std::vector< int > flybySequence,
                           const bool useTripTime = false );

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    std::string get_name( ) const;

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        if(useTripTime_ )
        {
            return 2u;
        }
        else
        {
            return 1u;
        }

    }

private:

    mgaproblem mgaObject_;

    const std::vector< std::vector< double > > problemBounds_;

    bool useTripTime_;
};

#endif // TUDAT_EXAMPLE_PAGMO_MULTIPLE_GRAVITY_ASSIST_H
