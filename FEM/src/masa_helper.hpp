#include "elements.hpp"
#include <./eigen3/Eigen/Dense> 
#include <masa.h>
#include <string>
#include <vector>

using namespace Eigen;

class Manufactured {
    private:
        Manufactured(){};
    public:
        Manufactured(const char* name, const char* analysis) ;
        double forcing(double x) ;
        double stiffness ;
        VectorXd exact_solution( std::vector<Node1d>& Nodes );
};
