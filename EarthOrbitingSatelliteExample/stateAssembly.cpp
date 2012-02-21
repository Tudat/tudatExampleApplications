/*!   Copyright (c) 2010-2012 Delft University of Technology.
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

#include "stateAssembly.h"

namespace earth_orbiting_satellite_example
{

//! Assemble state.
AssembledStateWithBodyIndices assembleState( ListOfStates listOfStates )
{
    // Declare assembled state based on the fact that each body state is a Cartesian state of
    // length 6.
    Eigen::VectorXd assembledState( listOfStates.size( ) * SIZE_OF_CARTESIAN_STATE );

    // Declare assembled state index.
    unsigned int assembledStateIndex = 0;

    // Declare associative map of bodies and their assembled state indices.
    BodyIndices assembledStateIndicesForBodies;

    // Assemble state.
    for ( ListOfStates::iterator iteratorStates = listOfStates.begin( );
          iteratorStates != listOfStates.end( ); iteratorStates++ )
    {
        assembledState.segment( assembledStateIndex, SIZE_OF_CARTESIAN_STATE )
                = iteratorStates->second;

        assembledStateIndicesForBodies[ assembledStateIndex ] = iteratorStates->first;

        assembledStateIndex += SIZE_OF_CARTESIAN_STATE;
    }

    // Return assembled state.
    return std::make_pair( assembledState, assembledStateIndicesForBodies );
}

//! Disassemble state.
ListOfStates disassembleState(
        AssembledStateWithBodyIndices assembledStateWithBodyIndices )
{
    // Declare size of Cartesian state.
    unsigned int SIZE_OF_CARTESIAN_STATE = 6;

    // Declare list of states.
    ListOfStates listOfStates;

    // Disassemble state.
    for ( BodyIndices::iterator iteratorBodyIndices
          = assembledStateWithBodyIndices.second.begin( );
          iteratorBodyIndices != assembledStateWithBodyIndices.second.end( );
          iteratorBodyIndices++ )
    {
        listOfStates[ iteratorBodyIndices->second ] =
                assembledStateWithBodyIndices.first.segment(
                                    iteratorBodyIndices->first, SIZE_OF_CARTESIAN_STATE );
    }

    // Return list of states.
    return listOfStates;
}

} // namespace earth_orbiting_satellite_example