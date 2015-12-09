#include "elements.hpp"
#include "solvers.hpp"
#include "masa_helper.hpp"
#include <./eigen3/Eigen/Dense>
#include <masa.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace MASA ;

int main() {

	int order_init = 2 ; 
	int nx_init = 256; 
	int quad_pts_init = 4 ; 
	double Xmin_init = 0 ; 
	double Xmax_init = 4 ; 
	double dirichlet_init = 0.0 ;
	int method_init = 2 ;
	double tol_init = 1e-14 ;
	int max_iter_init = 1000000 ;
    int report_interval = 1000 ;
    bool report = true ;
    bool with_masa = true ;
    std::function<double(double)> stiffness ;
    std::function<double(double)> forcing ;
    
    Manufactured mfg("Chris Test 2","heateq_1d_steady_const");
    
    if(with_masa){
        stiffness = [&mfg] (double x)->double{return mfg.stiffness;} ;
        forcing = [&mfg] (double x)->double{return mfg.forcing(x);} ;
    } else {
        stiffness = [] (double x)->double{return 1;} ;
        forcing = [] (double x)->double{return 1;} ;
    }
    
    Domain1d Domain(order_init, nx_init, quad_pts_init, Xmin_init, Xmax_init, dirichlet_init ) ;
	Domain.forcing = forcing ;
	Domain.stiffness = stiffness ; 
    Domain.with_masa = with_masa ;
    Domain.build_elements( ) ;
    if (order_init == 2) {
        Domain.addNodes(1);
    }

	Solver GaussSeidel(method_init, tol_init, max_iter_init, report_interval, report) ;

 	GaussSeidel.solution_init(Domain.Nodes.size()) ;

	Domain.build_Ab(&GaussSeidel.A, &GaussSeidel.b);
	Domain.add_constraints(GaussSeidel.A, GaussSeidel.b);

 	GaussSeidel.solve() ;
    VectorXd exact = mfg.exact_solution(Domain.Nodes) ;
   
    std::cout << (exact-GaussSeidel.x).norm() << std::endl ;

    std::ofstream f ;
    f.open("masa_validate_2nd_256.txt");
    f << "x,exact,approx" << std::endl ;
    for(int i = 0; i < Domain.Nodes.size(); i++){
        f << Domain.Nodes[i].coords << "," << exact(i) << ',' << GaussSeidel.x(i) << std::endl ;
    }
    f.close();
}
