#include <pagmo/pagmo.hpp>

using namespace pagmo;

struct problem_basic {
    // Mandatory, computes ... well ... the fitness
    vector_double fitness(const vector_double &x) const
    {
        return {x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]};
    }

    // Mandatory, returns the box-bounds
    std::pair<vector_double, vector_double> get_bounds() const
    {
        return {{-10, -10, -10, -10}, {10, 10, 10, 10}};
    }

};



int main()
{
    // 1 - Instantiate a pagmo problem constructing it from a UDP
    // (user defined problem).
    problem prob{problem_basic()};

    // 2 - Instantiate a pagmo algorithm
    algorithm algo{sade(100)};

    // 3 - Instantiate an archipelago with 16 islands having each 20 individuals
    archipelago archi{16, algo, prob, 20};

    // 4 - Run the evolution in parallel on the 16 separate islands 10 times.
    archi.evolve(10);

    // 5 - Wait for the evolutions to be finished
    archi.wait();

    // 6 - Print the fitness of the best solution in each island
    for (const auto &isl : archi) {
        print(isl.get_population().champion_x(), "\n");
    }
}

