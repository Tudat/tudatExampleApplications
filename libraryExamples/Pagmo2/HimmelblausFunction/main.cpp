#include <iostream>
#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/pso.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/sga.hpp"
#include "pagmo/island.hpp"
#include "pagmo/problem.hpp"
#include "himmelblaus.hpp"

int main( )
{      

    pagmo::random_device::set_seed( 12345 );

    pagmo::problem prob{ my_problem( -5, 5, -5, 5) };

    pagmo::algorithm algo{ pagmo::de( ) };

    pagmo::island isl = pagmo::island{ algo, prob, 100 };

    for( int i = 1; i <= 1000; i++ )
    {

        isl.evolve( );

        //isl.get_population( ).get_f()

        //std::cout << "Best x: " << isl.get_population().champion_x()[0] << std::endl;
        //std::cout << "Best y: " << isl.get_population().champion_x()[1] << std::endl;
        std::cout << "Minimum: " << std::setprecision( 16 ) << isl.get_population().champion_f()[0] << std::endl;
        while( isl.status()!=pagmo::evolve_status::idle )
            isl.wait();

    }

    std::cout << "Best x: " << isl.get_population().champion_x()[0] << std::endl;
    std::cout << "Best y: " << isl.get_population().champion_x()[1] << std::endl;
    std::cout << "Minimum: " << isl.get_population().champion_f()[0] << std::endl;

    return 0;

}
