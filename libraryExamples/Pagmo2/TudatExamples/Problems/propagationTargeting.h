/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_EXAMPLE_PAGMO_PROBLEM_PROPAGATION_TARGETING_HPP
#define TUDAT_EXAMPLE_PAGMO_PROBLEM_PROPAGATION_TARGETING_HPP

#include <utility>
#include <vector>


// Define the problem PaGMO-style
struct PropagationTargetingProblem {

    // Empty constructor
    PropagationTargetingProblem( ){ }

    PropagationTargetingProblem( const double altitudeOfPerigee, const double altitudeOfApogee,
                                 const double altitudeOfTarget, const double longitudeOfTarget );

    // Fitness: takes the value of the RAAN and returns the value of the closest distance from target
    std::vector<double> fitness(const std::vector<double> &x) const;

    // Boundaries of the problem set between 0 and (360) degrees
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const;

private:

    double altitudeOfPerigee_;
    double altitudeOfApogee_;
    double altitudeOfTarget_;
    double longitudeOfTarget_;

};

#endif // TUDAT_EXAMPLE_PAGMO_PROBLEM_PROPAGATION_TARGETING_HPP
