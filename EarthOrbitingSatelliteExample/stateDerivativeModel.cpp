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

#include "stateAssembly.h"
#include "stateDerivativeModel.h"

namespace earth_orbiting_satellite_example
{

//! Compute state derivative.
Eigen::VectorXd StateDerivativeModel::computeStateDerivative(
        const double independentVariable, const Eigen::VectorXd& stateVector )
{
    // Declare assembled state derivative of the same size as state vector.
    Eigen::VectorXd assembledStateDerivative = Eigen::VectorXd::Zero( stateVector.rows( ) );

    // Set loop counter.
    unsigned int loopCounter = 0;

    // Generate assembled state-derivative vector.
    for ( ListOfForces::iterator iteratorForForces = listOfForces_.begin( );
          iteratorForForces != listOfForces_.end( ); iteratorForForces++ )
    {
        assembledStateDerivative.segment( loopCounter, 3 )
                = stateVector.segment( loopCounter + 3, 3 );

        for ( unsigned int i = 0; i < iteratorForForces->second.size( ); i++ )
        {
            tudat::State temporaryState;
            temporaryState.state = stateVector.segment( loopCounter, 6 );

            // Compute force.
            iteratorForForces->second.at( i )->computeForce(
                        &temporaryState, independentVariable );

            // Compute force.
            assembledStateDerivative.segment( loopCounter + 3, 3 )
                    += iteratorForForces->second.at( i )->getForce( )
                    / iteratorForForces->first->getMass( );
        }

        loopCounter += SIZE_OF_CARTESIAN_STATE;
    }

    // Return assembled state derivative.
    return assembledStateDerivative;
}

} // namespace earth_orbiting_satellite_example
