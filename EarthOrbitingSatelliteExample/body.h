/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      121030    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef EARTH_ORBITING_SATELLITE_EXAMPLE_BODY_H
#define EARTH_ORBITING_SATELLITE_EXAMPLE_BODY_H

#include <map>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>

namespace satellite_example
{

//! Typedef for Eigen::Vector6d.
typedef Eigen::Matrix< double, 6, 1 > Vector6d;

//! Typedef for propagation history container.
typedef std::map< double, Eigen::VectorXd > PropagationHistory;

//! Test body class.
/*!
 * This class serves as an example of how a container can be constructed that stored state and
 * time information, which can be used in conjunction with acceleration models, the Cartesian state
 * derivative model, and the composite state derivative model. It should be noted that this class
 * should NOT be used as is in the rest of the library, and is only used in conjunction with unit
 * tests. Classes of this nature should go in user applications, since they are typically
 * application-specific, hence not general enough to be added to the Tudat libraries.
 * \tparam SpatialDimensions Spatial dimensions of position and velocity vectors. The state vector
 *          is taken to be twice the number of spatial dimensions. All vectors are represented by
 *          Eigen::Matrix-types.
 * \tparam DataType The data type of all of the internally stored state-related vectors and time.
 */
class Body
{
public:

    //! Constructor taking a state and a time.
    /*!
     * Constructor taking an input state and time. The input state is used internally to
     * set the current position (taken as a segment of the input state given by the indices
     * (0, SpatialDimensions)) and the current velocity (taken as a segment of the input state
     * given by the indices (SpatialDimensions, SpatialDimensions).
     * \param aState An input state vector.
     * \param aTime An input time.
     */
    Body( const Vector6d& aState, const double aTime )
        : currentState( aState ),
          currentPosition( aState.segment( 0, 3 ) ),
          currentVelocity( aState.segment( 3, 3 ) ),
          currentTime( aTime )
    { }

    //! Set current time and state.
    /*!
     * Sets the current time, position and current velocity internally based on the input
     * arguments. The current position is taken as a segment of the input state given by the
     * indices (0, SpatialDimensions)), and the current velocity is taken as a segment of the input
     * state given by the indices (SpatialDimensions, SpatialDimensions).
     * \param aTime An input time.
     * \param aState An input state vector.
     */
    void setCurrentTimeAndState( const double aTime, const Vector6d& aState )
    {
        currentTime = aTime;
        currentState = aState;
        currentPosition = aState.segment( 0, 3 );
        currentVelocity = aState.segment( 3, 3 );
    }

    //! Get current state.
    /*!
     * Returns the internally stored current state vector.
     * \return Current state.
     */
    Vector6d getCurrentState( ) { return currentState; }

    //! Get current position.
    /*!
     * Returns the internally stored current position vector.
     * \return Current position.
     */
    Eigen::Vector3d getCurrentPosition( ) { return currentPosition; }

    //! Get current velocity.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity.
     */
    Eigen::Vector3d getCurrentVelocity( ) { return currentVelocity; }

    //! Get current time.
    /*!
     * Returns the internally stored current time.
     * \return Current time.
     */
    double getCurrentTime( ) { return currentTime; }

protected:

private:

    //! Current state.
    Vector6d currentState;

    //! Current position.
    Eigen::Vector3d currentPosition;

    //! Current position.
    Eigen::Vector3d currentVelocity;

    //! Current time.
    double currentTime;
};

//! Typedef for shared-pointer to body.
typedef boost::shared_ptr< Body > BodyPointer;

//! Typedef for map of satellites with associated list of acceleration models.
typedef std::map< BodyPointer, std::vector< tudat::basic_astrodynamics::
AccelerationModel3dPointer > > ListOfSatellites;

//! Data updater class.
/*!
 * This class is used to update the states of a list of satellites. It serves as a data repository
 * to communicate the state generated by the numerical integrator with the acceleration models
 * used. It works with composite states.
 */
class DataUpdater
{
public:

    //! Constructor taking list of satellites.
    /*!
     * This constructor takes a list of satellites, for which the states are managed and updated by
     * the class.
     * \param someSatellites List of satellites.
     */
    DataUpdater( const ListOfSatellites someSatellites )
        : satellites( someSatellites )
    { }

    //! Update body data.
    /*!
     * Updates body data by decomposing the composite state into its constituent states, and
     * storing the results in the Body objects associated with each satellite.
     * \param time Current time.
     * \param compositeState Composite state (Eigen::VectorXd-type).
     */
    void updateBodyData( const double time, const Eigen::VectorXd compositeState )
    {
        unsigned int counter = 0;

        for ( ListOfSatellites::const_iterator it = satellites.begin( );
              it != satellites.end( ); it++ )
        {
            it->first->setCurrentTimeAndState( time, compositeState.segment( counter * 6, 6 ) );
            counter++;
        }
    }

protected:

private:

    //! List of satellites.
    /*!
     * List of satellites provided to the class, containing the body objects for each satellite,
     * and the associated lists of acceleration models.
     */
    const ListOfSatellites satellites;
};

} // namespace example_satellite

#endif // EARTH_ORBITING_SATELLITE_EXAMPLE_BODY_H
