/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "multipleGravityAssist.h"

MultipleGravityAssist::MultipleGravityAssist( std::vector< std::vector< double > > &bounds,
                                              const bool useTripTime ) :
    problemBounds_( bounds ), useTripTime_( useTripTime ) {

    std::vector< int > flybySequence = { 3, 4 };
    mgaObject_.type = total_DV_rndv;
    mgaObject_.sequence = flybySequence;
    mgaObject_.rev_flag.resize( flybySequence.size( ) );

    for( unsigned int i = 0; i < flybySequence.size( ); i++ )
    {
        mgaObject_.rev_flag[ i ] = 0;
    }

    mgaObject_.Isp = 300.0;
    mgaObject_.mass = 8000.0;
    mgaObject_.DVlaunch = 0.0;
}


//! Descriptive name of the problem
std::string MultipleGravityAssist::get_name() const {
    return "MGA transfer trajectory";
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > MultipleGravityAssist::get_bounds() const {

    return { problemBounds_[0], problemBounds_[1] };
}

//! Implementation of the fitness function (return delta-v)
std::vector<double> MultipleGravityAssist::fitness( const std::vector<double> &xv ) const{

            std::vector<double> rp;
            std::vector<double> DV;
            double obj_funct;
            int result = MGA( xv, mgaObject_, rp, DV, obj_funct );
            return { obj_funct };
}


