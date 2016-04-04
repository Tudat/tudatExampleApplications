/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H
#define TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H

#include <vector>

#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

#include <pagmo/src/config.h>
#include <pagmo/src/serialization.h>
#include <pagmo/src/types.h>
#include <pagmo/src/problem/base.h>

#include <Eigen/Core>

namespace pagmo
{

namespace problem
{

typedef Eigen::Matrix< double, 6, 1 > StateType;

//! Test function for a new interplanetary trajectory class in Tudat
class __PAGMO_VISIBLE EarthMarsTransfer : public base
{
  public:
    EarthMarsTransfer( const std::vector< std::vector< double > > problemBounds );
    base_ptr clone( ) const;
    std::string get_name( ) const;
    const std::vector< std::vector< double > > problemBounds_;

  protected:
    void objfun_impl( fitness_vector& f, const decision_vector& xv ) const;

  private:
    StateType getPlanetPosition( const double date, const std::string planetName ) const;

    friend class boost::serialization::access;
    template< class Ar > void serialize( Ar &ar, const unsigned int )
    {
        ar & boost::serialization::base_object< base >( *this );
    }
};

} // namespace problem;

} // namespace pagmo

// BOOST_SERIALIZATION_ASSUME_ABSTRACT( pagmo::problem::EarthMarsTransfer );
namespace boost {

// Serialization for Eigen vectors and matrices
template<class Ar, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void serialize(
    Ar &ar,
    Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,
    const unsigned int file_version
)
{
    for( size_t i = 0; static_cast< unsigned int >( i )< t.size( ); i++ )
        ar & t.data( )[ i ];
}

// Because we have a non-empty constructor we need to be careful to
// specify the save and load procedure for serialization.
namespace serialization
{

template<class Ar>
inline void save_construct_data( Ar &ar, const pagmo::problem::EarthMarsTransfer * t,
    const unsigned int file_version ){

    // save data required to construct instance
    ar << t->problemBounds_;
}

template<class Ar>
inline void load_construct_data( Ar &ar, pagmo::problem::EarthMarsTransfer * t,
    const unsigned int file_version ){

    // retrieve data from archive required to construct new instance
    std::vector< std::vector< double > > problemBounds;
    ar >> problemBounds;
    // invoke inplace constructor to initialize instance of my_class
    ::new( t )pagmo::problem::EarthMarsTransfer( problemBounds );
}

} // namespace serialization;

} // namespace boost



BOOST_CLASS_EXPORT_KEY( pagmo::problem::EarthMarsTransfer );

#endif // TUDAT_EXAMPLE_PAGMO_PROBLEM_EARTH_MARS_TRANSFER_H
