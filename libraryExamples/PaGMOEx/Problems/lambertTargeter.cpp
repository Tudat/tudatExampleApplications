/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "PaGMOEx/Problems/lambertTargeter.h"


namespace pagmo { namespace problem {

EarthMarsTransfer::EarthMarsTransfer( const std::vector< std::vector< double > > problemBounds ) :
    base( problemBounds[0], problemBounds[1], 0, 2 );

//! Clone method.
base_ptr EarthMarsTransfer::clone( ) const {
        return base_ptr(new generalLearning(*this));
}

//! Implementation of the objective function.
void generalLearning::objfun_impl(fitness_vector &f, const decision_vector &xv) const{


  StateType EarthMarsTransfer::getPlanetPosition( const double date,
                                                  const std::string planetName ) const {
    using tudat::orbital_element_conversions::convertKeplerianToModifiedEquinoctialElements;
    using tudat::orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly;
    StateType stateKepl, stateCart;
    double n  = 1.991e-07;
    if( planetName == "Earth" ){
        stateKepl << 1.4960e+11, 1.6717e-02, 0.0, 5.0198e+00, 3.0614e+00, 4.4961e+00, 5.4000e+04;
    else
        stateKepl << 2.2794e+11, 9.3412e-02, 0.0, 5.8650e+00, 8.6531e-01, 5.7567e+00, 5.1544e+04;
    }
    ConvertMeanAnomalyToEccentricAnomaly meanToTrueConverter( stateKepl( 1 ),
        fmod( stateKepl( 5 ) + ( date - stateKepl( 6 ) ) * 86400. * n, 2.*M_PI ) );
    stateKepl( 5 ) = meanToTrueConverter.convert();
    stateCart = convertKeplerianToCartesianCoordinates( stateKepl , false );
    return stateCart;
}

}} // namespace problem; namespace pagmo
