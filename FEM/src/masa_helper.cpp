#include "masa_helper.hpp"
#include "elements.hpp"
#include <./eigen3/Eigen/Dense>
#include <masa.h>
#include <string>
#include <vector>
#include <iostream>

using namespace Eigen;

Manufactured::Manufactured(const char* name, const char* analysis) {
    
    masa_init ( "chris" , "heateq_1d_steady_const" );
    masa_init_param();
    stiffness = masa_get_param("k_0");
}

double Manufactured::forcing(double x) {
    double out;
    out = masa_eval_1d_source_t(x) ;
    return out ;
}

VectorXd Manufactured::exact_solution( std::vector<Node1d>& Nodes ) {
    VectorXd exact = VectorXd::Zero(Nodes.size()) ;

    for (int i = 0; i < Nodes.size(); i++) {
        exact(i) = masa_eval_1d_exact_t( Nodes[i].coords ) ;
    }
    return exact ;
}




