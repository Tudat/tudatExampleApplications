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
