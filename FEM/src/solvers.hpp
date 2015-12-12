#ifndef SOLVERS_H
#define SOLVERS_H
#include "elements.hpp"
#include "./eigen3/Eigen/Dense"
#include <iostream>
#include <fstream>

class Solver {
public:
    std::string file_out ;
    int iter, report_interval;
	double tol ;
	int max_iter, method, n;
	bool report, with_masa;
	double res ;
	Eigen::MatrixXd A ;
	Eigen::VectorXd b ; 
	Eigen::VectorXd x ;

	Solver(){} ;
	Solver(int method_init, double tol_init, int max_iter_init, int report_interval_init, bool report_init);
	void solution_init(int n_eq) ;
	void solve() ;
	void jacobi() ;
	void gaussSeidel() ;
    void Householder() ;
    void CG() ;
    void output(Domain1d &Domain);
 
} ;

#endif
