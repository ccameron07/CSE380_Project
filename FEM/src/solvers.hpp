#ifndef SOLVERS_H
#define SOLVERS_H
#include "./eigen3/Eigen/Dense"
#include <iostream>

class Solver {
private:

	Solver(){} ;
	double tol ;
	int max_iter, method, n;

public:
	int iter ;
	double res ;
	Eigen::MatrixXd A ;
	Eigen::VectorXd b ; 
	Eigen::VectorXd x ;

	Solver(int method_init, double tol_init, int max_iter_init) ;
	void solution_init(int n_eq) ;
	void solve() ;
	void jacobi() ;
	void gaussSeidel() ;
 
} ;

#endif
