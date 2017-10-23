#include<iostream>
#include"pagmo/algorithms/sade.hpp"
#include"pagmo/island.hpp"
#include"pagmo/problem.hpp"
#include"himmelblaus.hpp"

int main( )
{      

    pagmo::random_device::set_seed( 12345 );

    pagmo::problem prob{ my_problem( 0, 5, 0, 5) };

    pagmo::algorithm algo{ pagmo::sade( ) };

    pagmo::island isl = pagmo::island{ algo, prob, 128 };

    for( int i = 1; i <= 10000; i++ )
    {

        isl.evolve( );

        while( isl.status()!=pagmo::evolve_status::idle )
            isl.wait();

    }

    std::cout << "Best x: " << isl.get_population().champion_x()[0] << std::endl;
    std::cout << "Best y: " << isl.get_population().champion_x()[1] << std::endl;
    std::cout << "Minimum: " << isl.get_population().champion_f()[0] << std::endl;

    return 0;

}
