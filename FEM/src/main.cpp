#include <vector>
#include "./eigen3/Eigen/Dense"
#include "elements.hpp"
#include "solvers.hpp"
#include <iostream>

int main() {

	int order_init = 1 ; 
	int nx_init = 100 ; 
	int quad_pts_init = 3 ; 
	double Xmin_init = 0.0 ; 
	double Xmax_init = 3.0 ; 
	double dirichlet_init = 0.0 ;
	int method_init = 1 ;
	double tol_init = 1e-11 ;
	int max_iter_init = 10000 ;

	Domain1d Domain(order_init, nx_init, quad_pts_init, Xmin_init, Xmax_init, dirichlet_init ) ;
	Domain.build_elements( ) ;

	Solver GaussSeidel(method_init, tol_init, max_iter_init) ;

 	GaussSeidel.solution_init(Domain.Nodes.size()) ;

	Domain.build_Ab(&GaussSeidel.A, &GaussSeidel.b);
	Domain.add_constraints(GaussSeidel.A, GaussSeidel.b);

 	GaussSeidel.solve() ;
 	
}