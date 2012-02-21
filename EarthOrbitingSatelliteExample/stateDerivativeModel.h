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
 *      120220    K. Kumar          File created.
 *      120221    K. Kumar          Rewritten for new application architecture.
 *
 *    References
 *
 */

#ifndef STATE_DERIVATIVE_MODEL_H
#define STATE_DERIVATIVE_MODEL_H

#include <Eigen/Core>
#include <map>
#include <vector>
#include <Tudat/Astrodynamics/BasicAstrodynamics/forceModel.h>
#include <Tudat/Astrodynamics/Bodies/body.h>

namespace earth_orbiting_satellite_example
{

//! State derivative class.
/*!
 * Definition of the state derivative class for an N-body system. This Cartesian state derivative
 * (derivative of [x,y,z,xdot,ydot,zdot] computes the orbiting motion of a system of bodies, given
 * a set of forces acting on them.
 */
class StateDerivativeModel
{
public:

    //! Typedef for list of forces per body.
    typedef std::map< tudat::Body*, std::vector< tudat::ForceModel* > > ListOfForces;

    //! Default constructor.
    /*!
     * Default constructor, requiring list of forces per body.
     * \param listOfForces Associative container of list of forces per body.
     */
    StateDerivativeModel( ListOfForces listOfForces ) : listOfForces_( listOfForces ) { }

    //! Compute state derivative.
    /*!
     * Computes the state derivative given to the integrator. This is an implementation of the
     * function-type required by the NumericalIntegrator base class in Tudat Core.
     * \param independentVariable Independant variable.
     * \param stateVector State vector.
     * \return State-derivative vector.
     */
    Eigen::VectorXd computeStateDerivative( const double independentVariable,
                                            const Eigen::VectorXd& stateVector );

protected:

private:

    //! List of forces per body.
    /*!
     * List of forces per body.
     */
    ListOfForces listOfForces_;
};

} // namespace earth_orbiting_satellite_example

#endif // STATE_DERIVATIVE_MODEL_H