#include "elements.hpp"
#include "solvers.hpp"
#include "masa_helper.hpp"
#include <./eigen3/Eigen/Dense>
#include <masa.h>
#include <iostream>
#include <vector>

using namespace MASA ;

int main() {

	int order_init = 1 ; 
	int nx_init = 100 ; 
	int quad_pts_init = 3 ; 
	double Xmin_init = -1.121997376282069 ; 
	double Xmax_init =  1.121997376282069 ; 
	double dirichlet_init = 0.0 ;
	int method_init = 1 ;
	double tol_init = 1e-14 ;
	int max_iter_init = 100000 ;
    int report_interval = 1000 ;
    bool report = true ;
    
    masa_init("Chris' Test","heateq_1d_steady_const");
    masa_display_param();

    Manufactured mfg("Chris Test 2","heateq_1d_steady_const");
    std::cout << mfg.stiffness << std::endl;
    std::cout << mfg.forcing(0.5) << std::endl;
    //std::function<double(double)> stiffness1 = [&mfg] (double x)->double{return mfg.forcing(x);};

    std::function<double(double)> stiffness = [&mfg] (double x)->double{return mfg.stiffness;} ;
    std::function<double(double)> forcing = [&mfg] (double x)->double{return mfg.forcing(x);} ;
	Domain1d Domain(order_init, nx_init, quad_pts_init, Xmin_init, Xmax_init, dirichlet_init ) ;
	Domain.forcing = forcing ;
	Domain.stiffness = stiffness ; 
    Domain.build_elements( ) ;

	Solver GaussSeidel(method_init, tol_init, max_iter_init, report_interval, report) ;

 	GaussSeidel.solution_init(Domain.Nodes.size()) ;

	Domain.build_Ab(&GaussSeidel.A, &GaussSeidel.b);
	Domain.add_constraints(GaussSeidel.A, GaussSeidel.b);

 	GaussSeidel.solve() ;
    std::cout << GaussSeidel.x << std::endl;
    std::cout << "-------------------" << std::endl;
    mfg.exact_solution(Domain.Nodes) ;
    for(int i = 0; i < Domain.Nodes.size(); i++){
        std::cout << Domain.Nodes[i].coords << std::endl ;
    }
}
