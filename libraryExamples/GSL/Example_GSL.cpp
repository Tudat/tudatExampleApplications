/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include <iostream> // cout sometimes needs this

#include <cstdio>
#include <vector>
#include <fstream>
#include <iomanip> // setw

#include <Eigen/Core>

// GSL
#include <gsl/gsl_integration.h>    // numerical quadrature
#include <gsl/gsl_qrng.h>           // GSL Quasi Random generators -> Sobol

//// ============================ Header / Declarations ====================== ////

using namespace std;

// function to be integrated
double my_integrand (double x, void *params)
{
  // Mathematica form:  Log[alpha*x]/Sqrt[x]

  // The next line recovers alpha from the passed params pointer
  double alpha = *(double *) params;

  return (log (alpha * x) / sqrt (x));
}

//// =================== GSL Numerical quadrature ==================== ////
/// https://www.physics.ohio-state.edu/~ntg/780/gsl_examples/qags_test.cpp
/// https://www.physics.ohio-state.edu/~ntg/780/gsl_examples/
void gsl_integrator_example(){
    gsl_integration_workspace *work_ptr
      = gsl_integration_workspace_alloc (1000);

    double lower_limit = 0;	/* lower limit a */
    double upper_limit = 1;	/* upper limit b */
    double abs_error = 1.0e-8;	/* to avoid round-off problems */
    double rel_error = 1.0e-8;	/* the result will usually be much better */
    double result;		/* the result from the integration */
    double error;			/* the estimated error from the integration */

    double alpha = 1.0;		// parameter in integrand
    double expected = -4.0;	// exact answer

    gsl_function My_function;
    void *params_ptr = &alpha;

    My_function.function = &my_integrand;
    My_function.params = params_ptr;

    gsl_integration_qags (&My_function, lower_limit, upper_limit,
              abs_error, rel_error, 1000, work_ptr, &result,
              &error);

    cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
    cout.precision (18);		// 18 digits in doubles

    int width = 20;  // setw width for output
    cout << "result          = " << setw(width) << result << endl;
    cout << "exact result    = " << setw(width) << expected << endl;
    cout << "estimated error = " << setw(width) << error << endl;
    cout << "actual error    = " << setw(width) << result - expected << endl;
    cout << "intervals =  " << work_ptr->size << endl;
}

//// =================== GSL Sobol Sampler ==================== ////
std::vector< Eigen::VectorXd > GSL_Sobol_Sampler(const int Dimension, const int N){
    std::vector< Eigen::VectorXd > SobolX(N);

    Eigen::VectorXd x( Dimension );
    double v[ Dimension ];

    gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, Dimension );

    for(int j = 0 ; j < N ; j++){ // Loop over samples
        gsl_qrng_get (q, v); // Generate sobol [0,1]

        for(int i = 0 ; i < Dimension ; i++){ // Fill vector
            x(i) = v[i] ;
        }
        SobolX[j] = x ; // Save vector
    }

    gsl_qrng_free (q); // don't know why?

    return SobolX;
}


//// =================== Start Main code ==================== ////
int main (void)
{
    std::cout << "========= EXAMPLES USING GSL LIBRARY =========" << std::endl << std::endl ;
    std::cout << "Numerical quadrature example using GSL" << std::endl;

    gsl_integrator_example();

    std::cout << std::endl;
    std::cout << "Sobol sampling example using GSL" << std::endl;

    std::vector< Eigen::VectorXd > Samples = GSL_Sobol_Sampler(3,10) ;

    std::cout << "Samples:" << std::endl;
    for(unsigned int i = 0 ; i < Samples.size() ; i++){
        std::cout << Samples[i] << std::endl << std::endl;
    }

  return 0;
}
