/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef PAGMO_PROBLEM_SAMPLE_RETURN_H
#define PAGMO_PROBLEM_SAMPLE_RETURN_H


#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

#include <pagmo/src/config.h>
#include <pagmo/src/serialization.h>
#include <pagmo/src/types.h>
#include <pagmo/src/problem/base.h>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h>


namespace pagmo { namespace problem {

typedef Eigen::Matrix< double, 6, 1 > StateType;

//! Test function for a new interplanetary trajectory class in Tudat
class __PAGMO_VISIBLE sampleReturn : public base
{
  public:
    sampleReturn( const StateType);
    base_ptr clone() const;

  protected:
    void objfun_impl(fitness_vector &, const decision_vector &) const;

  private:
    StateType getEarthPosition( double jd );
    StateType getMarsPosition( double jd );
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & boost::serialization::base_object<base>(*this);
    }
}
}} // namespace problem; namespace pagmo

// BOOST_SERIALIZATION_ASSUME_ABSTRACT( pagmo::problem::sampleReturn );
namespace boost {

// Serialization for Eigen vectors and matrices
template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void serialize(
    Archive & ar,
    Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,
    const unsigned int file_version
)
{
    for(size_t i=0; static_cast< unsigned int >(i)<t.size(); i++)
        ar & t.data()[i];
}

// Because we have a non-empty constructor we need to be careful to
// specify the save and load procedure for serialization.
namespace serialization {

template<class Archive>
inline void save_construct_data(
    Archive & ar, const pagmo::problem::neuroTransfer * t,
    const unsigned int file_version ){

    // save data required to construct instance
    ar << t->initialState << t->finalState << t->weights << t->constraints
       << t->bounds << t->dims << t->params
       << t->ptype;
}

template<class Archive>
inline void load_construct_data( Archive & ar,
    pagmo::problem::sampleReturn * t,
    const unsigned int file_version ){

    // retrieve data from archive required to construct new instance
    pagmo::problem::StateType initialState;
    pagmo::problem::StateType finalState;
    pagmo::problem::StateType weights;
    std::vector< std::vector< int > > constraints;
    std::vector< double > bounds;
    std::vector< int > dims;
    std::map< std::string, double > params;
    pagmo::problem::neuroTransfer::ProblemType ptype;
    ar >> initialState >> finalState >> weights >> constraints
       >> bounds >> dims >> params
       >> ptype;
    // invoke inplace constructor to initialize instance of my_class
    ::new(t)pagmo::problem::sampleReturn(
        initialState, finalState, weights, constraints,
        bounds, dims, params, ptype );
}
}
}
