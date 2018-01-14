/*    Copyright (c) 2010-2018, Delft University of Technology
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
 *      101111    K. Kumar          File created.
 *      110113    K. Kumar          Scenario updated to use latest version of code; added file
 *                                  header and footer.
 *      110202    K. Kumar          Scenario updated to use latest version of code.
 *      110216    K. Kumar          Migrated to applications namespace.
 *      110217    K. Kumar          Function name changed.
 *      110815    K. Kumar          Updated with mass of Asterix.
 *      111024    K. Kumar          Modified to be executable program with main-function as
 *                                  suggested by M. Persson.
 *      120221    K. Kumar          Rewrote application from scratch; now propagates two
 *                                  satellites.
 *      120502    K. Kumar          Updated code to use shared pointers.
 *      121030    K. Kumar          Updated code to use new state derivative models.
 *      130107    K. Kumar          Updated license in file header.
 *      130225    K. Kumar          Updated gravitational acceleration model references; renamed
 *                                  file; fixed error in assigning Obelix state derivative model;
 *                                  made variables const-correct.
 *
 *      160322    R. Hoogendoorn    GSL Library example
 *
 *    References
 *
 *    https://github.com/ampl/gsl
 *    http://www.gnu.org/software/gsl/
 *
 *    Notes
 *
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

#include <Eigen/Core>

#include <gsl/gsl_integration.h>    // numerical quadrature
#include <gsl/gsl_qrng.h>           // GSL Quasi Random generators -> Sobol

//! Function to be integrated
double my_integrand( double independentVariable, void *parameters )
{
    // The next line recovers alpha from the passed parameters pointer

    double alpha = *(double *)parameters;

    return ( std::log ( alpha * independentVariable ) / std::sqrt ( independentVariable ) );
}

//// =================== GSL Numerical quadrature ==================== ////
/// https://www.physics.ohio-state.edu/~ntg/780/gsl_examples/qags_test.cpp
/// https://www.physics.ohio-state.edu/~ntg/780/gsl_examples/
void gsl_integrator_example( )
{
    gsl_integration_workspace* work_ptr = gsl_integration_workspace_alloc( 1000 );

    double result;		/* the result from the integration */
    double error;			/* the estimated error from the integration */

    // Define limits and tolerances
    double lower_limit = 0.0;	/* lower limit a */
    double upper_limit = 1.0;	/* upper limit b */
    double abs_error = 1.0E-8;	/* to avoid round-off problems */
    double rel_error = 1.0E-8;	/* the result will usually be much better */

    // Parameter in integrand
    double alpha = 1.0;

    // Exact answer
    double expected = -4.0;

    // Define integration function
    gsl_function my_function;
    void *parameters_ptr = &alpha;
    my_function.function = &my_integrand;
    my_function.params = parameters_ptr;

    // PErform integration.
    gsl_integration_qags(
                &my_function, lower_limit, upper_limit,
                abs_error, rel_error, 1000, work_ptr, &result,
                &error );

    std::cout.setf ( std::ios::fixed, std::ios::floatfield );	// output in fixed format
    std::cout.precision ( 18 );		// 18 digits in doubles

    int width = 20;  // setw width for output
    std::cout << "result          = " << std::setw( width ) << result << std::endl;
    std::cout << "exact result    = " << std::setw( width ) << expected << std::endl;
    std::cout << "estimated error = " << std::setw( width ) << error << std::endl;
    std::cout << "actual error    = " << std::setw( width ) << result - expected << std::endl;
    std::cout << "intervals =  " << work_ptr->size << std::endl;
}

//// =================== GSL Sobol Sampler ==================== ////
std::vector< Eigen::VectorXd > gsl_Sobol_Sampler( const int dimension, const int numberOfSamples )
{
    // Initialize vector list
    std::vector< Eigen::VectorXd > sobolX( numberOfSamples );

    Eigen::VectorXd currentVector( dimension );
    double randomArray[ dimension ];

    // Create random number generator
    gsl_qrng* randomNumberGenerator  = gsl_qrng_alloc( gsl_qrng_sobol, dimension );

    // Loop over samples
    for( int j = 0 ; j < numberOfSamples ; j++ )
    {
        // Generate sobol [0,1]
        gsl_qrng_get( randomNumberGenerator, randomArray );

        // Fill Eigen vector with results from raw array.
        for( int i = 0 ; i < dimension ; i++ )
        {
            currentVector( i ) = randomArray[ i ] ;
        }

        // Save vector
        sobolX[ j ] = currentVector ;
    }

    // Deallocate variables
    gsl_qrng_free( randomNumberGenerator );
    return sobolX;
}


//// =================== Start Main code ==================== ////
int main( )
{
    std::cout << "========= EXAMPLES USING GSL LIBRARY =========" << std::endl << std::endl ;
    std::cout << "Numerical quadrature example using GSL" << std::endl;

    // Run integration example.
    gsl_integrator_example( );

    std::cout << std::endl;
    std::cout << "Sobol sampling example using GSL" << std::endl;

    // Run sobol sampler example.
    std::vector< Eigen::VectorXd > samples = gsl_Sobol_Sampler( 3, 10 ) ;

    std::cout << "Samples:" << std::endl;
    for( unsigned int i = 0; i < samples.size( ) ; i++)
    {
        std::cout << samples[ i ] << std::endl << std::endl;
    }

    return 0;
}
