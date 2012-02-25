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
 *      120221    K. Kumar          File created.
 *
 *    References
 *
 */

#ifndef STATE_ASSEMBLY_H
#define STATE_ASSEMBLY_H

#include <Eigen/Core>
#include <map>
#include <utility>
#include "stateDerivativeModel.h"

namespace earth_orbiting_satellite_example
{

//! Size of Cartesian state vector for a single body.
const static unsigned int SIZE_OF_CARTESIAN_STATE = 6;

//! Typedef for list of states, per body.
typedef std::map< tudat::Body*, Eigen::VectorXd > ListOfStates;

//! Typedef for list of body indices in assembled state.
typedef std::map< int, tudat::Body* > BodyIndices;

//! Typedef for assembled state and body indices.
typedef std::pair< Eigen::VectorXd, BodyIndices > AssembledStateWithBodyIndices;

//! Assemble state.
/*!
 * Assembles state, given a list of states per body. This function produces a single state vector
 * used in conjunction with a numerical integrator, from a collection of states, associated with
 * different bodies.
 * \param listOfStates List of states associative with different bodies.
 * \return Assembled state and body indices relating state of body to position in assembled state.
 */
AssembledStateWithBodyIndices assembleState( ListOfStates listOfStates );

//! Disassemble state.
/*!
 * Disassembles state, given an associative container with the assembled stated, and the body
 * indices that relates the position of the state of a body to the index in the assembled state.
 * \param assembledStateWithBodyIndices Assembled state and associated indices of state of bodies
 *          within assembled state.
 * \return List of states, associated with different bodies.
 */
ListOfStates disassembleState( AssembledStateWithBodyIndices assembledStateWithBodyIndices );

} // namespace earth_orbiting_satellite_example

#endif // STATE_ASSEMBLY_H
