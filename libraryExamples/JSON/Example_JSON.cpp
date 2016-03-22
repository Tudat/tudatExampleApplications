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
 *      160322    R. Hoogendoorn    JSON Library example
 *
 *    References
 *
 *    https://github.com/open-source-parsers/jsoncpp.git
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

// JSon
#include <fstream>
#include <json/value.h>
#include <json/json.h>
#include <string>

//// ============================ Header / Declarations ====================== ////

//// Read config file and return Json::Value member
Json::Value ReadConfigFile(std::string folder, std::string filename){
    // Construct Json Database
    Json::Value settings;

    // Read file and overwrite Json database
    std::string filepath = folder+filename ;
    try{
        std::ifstream config_doc(filepath.c_str(), std::ifstream::binary);
        config_doc >> settings; // write file to settings
        std::cout << "Config file : " << settings["Name"].asString() << " (" <<settings["Date"][0]
                  << "-" <<settings["Date"][1] << "-"<<settings["Date"][2] <<") Version: "
                  << settings["Version"].asString() << std::endl;
    }
    catch(int ErrorNo){
        std::cerr << "Cannot read config file or write to Json::Value. Error No " << ErrorNo << std::endl;
    }

    return settings;
}

//// Read Json Array and return as Eigen Matrix
Eigen::MatrixXd ReadJsonMatrix(Json::Value Array){
    int rows = Array.size();
    int cols = Array[0].size();

    Eigen::MatrixXd Matrix(rows,cols) ;

    for(int i = 0 ; i<rows ; i++){
        for(int j=0 ; j<cols ; j++){
            Matrix(i,j) = Array[i][j].asDouble() ;
        }
    }
    return Matrix;
}

//// Read Json Array and return as Eigen Vector
Eigen::VectorXd ReadJsonVector(Json::Value Array){
    int size = Array.size();

    Eigen::VectorXd Vector(size);
    for(int i = 0 ; i < size ; i++){
        Vector(i) = Array[i].asDouble();
    }
    return Vector;
}

//// Read Json Array and return as Std Vector
std::vector<double> ReadJsonStdVector(Json::Value Array){
    int size = Array.size();

    std::vector<double> Vector(size);
    for(unsigned int i = 0 ; i < size ; i++){
        Vector[i] = Array[i].asDouble();
    }
    return Vector;
}


//// =================== Start Main code ==================== ////
int main (void)
{
    std::cout << "========= EXAMPLES USING JSON LIBRARY =========" << std::endl << std::endl ;

    // Open config file
    // save path of cpp file
    std::string cppPath( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string folder = cppPath.substr( 0, cppPath.find_last_of("/\\")+1);

    std::string filename ("Example_JSON.json");
    Json::Value settings = ReadConfigFile(folder,filename);

    // Read and cout values from "Example_JSON.json"
    std::cout << "Matrix A: " << ReadJsonMatrix(settings["A"]) << std::endl << std::endl;
    std::cout << "Vector X0: " << ReadJsonVector(settings["X0"]) << std::endl << std::endl;

    std::cout << "value = " << settings["value"].asDouble() << "  as int: "
              << settings["value"].asInt() << std::endl;
    std::cout << "value = " << settings["intvalue"].asInt() << std::endl;

    if(settings["tf"].asBool() == true){
        std::cout << "tf = true " << std::endl;
    }
    else{
        std::cout << "tf = false " << std::endl;
    }

  return 0;
}
