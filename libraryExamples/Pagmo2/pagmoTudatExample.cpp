/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


/* A satellite is going to be launched on an elliptical orbit with min and max altitudes
respectively 180 km and 40000 km with inclination i=0 (and argument of perigee omega = 0).
A fixed target is set on the same plane at an altitude of 35000 km, fixed latitude of 30 deg.
What is the value of the RAAN for which the satellite achieves a minimum
approach distance from the target?*/

#include <pagmo/pagmo.hpp>
#include "Pagmo2/Problems/easy_problem.hpp"

using namespace pagmo;

int main()
{

    //Set seed for reproducible results
    pagmo::random_device::set_seed(255);

    //Load spice kernels
    std::string kernelsPath = tudat::input_output::getSpiceKernelPath( );
    tudat::spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    tudat::spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");

    // 1 - Instantiate a pagmo problem constructing it from a UDP
    // (user defined problem).
    problem prob{easy_problem()};

    // 2 - Instantiate a pagmo algorithm
    algorithm algo{sade(10)};

    // 3 - Instantiate an archipelago with 16 islands having each 20 individuals
    archipelago archi{1, algo, prob, 16};

    // 4 - Run the evolution in parallel on the 16 separate islands 10 times.
    archi.evolve(20);

    // 5 - Wait for the evolutions to be finished
    archi.wait();

    // 6 - Print the decision and the fitness of the best solution in each island
    for (const auto &isl : archi) {
        print(isl.get_population().champion_x(), "Selected RAAN (deg)\n");
        print(isl.get_population().champion_f(), "Minimum separation from target (meters)\n");
    }
}

