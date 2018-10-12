#include <vector>
#include <iostream>
#include <functional>
#include <memory>

#include <Tudat/Basics/basicTypedefs.h>

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>
double testFunction( const double input, const double input2 )
{
    return 2.0 * input + input2;
}

class TestClass
{
public:
    TestClass( const double a ):a_( a ){ }

    double a_;
};

inline Eigen::Vector6d concatenateForceAndMomentCoefficients(
        const std::function< Eigen::Vector3d( const std::vector< double >& ) >&
        forceCoefficientFunction,
        const std::function< Eigen::Vector3d( const std::vector< double >& ) >&
        momentCoefficientFunction,
        const std::vector< double >& independentVariables )
{
    return ( Eigen::Vector6d( ) << forceCoefficientFunction( independentVariables ),
             momentCoefficientFunction( independentVariables ) ).finished( );
}


int main( )
{
   std::function< double( const double ) > myFunction =
           std::bind( &testFunction, std::placeholders::_1, 3.0 );
   std::cout<<myFunction( 1.0 )<<std::endl;
   std::cout<<( myFunction == nullptr )<<std::endl;

   std::function< double( const double ) > myFunction2;
   std::cout<<( myFunction2 == nullptr )<<std::endl;

   std::shared_ptr< TestClass > test = std::make_shared< TestClass >( 4.0 );


   std::function< Eigen::Vector3d( const std::vector< double >& ) > forceCoefficientFunction;
   std::function< Eigen::Vector3d( const std::vector< double >& ) > momentCoefficientFunction;


   std::function< Eigen::Vector6d( const std::vector< double >& ) > coefficientFunction;

   coefficientFunction = std::bind(
              &concatenateForceAndMomentCoefficients, forceCoefficientFunction, momentCoefficientFunction,
               std::placeholders::_1 );


}
