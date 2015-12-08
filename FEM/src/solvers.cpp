#include "./eigen3/Eigen/Dense"
#include <iostream>
#include "solvers.hpp"

using namespace Eigen;

Solver::Solver(int method_init, double tol_init, int max_iter_init, int report_interval_init, bool report_init) {
  method = method_init ;
  tol = tol_init ;
  max_iter = max_iter_init ;
  report_interval = report_interval_init ;
  report = report_init ;
}
    
void Solver::gaussSeidel(){

    VectorXd x1 = VectorXd::Zero(n) ;
    VectorXd C(n) ;
    MatrixXd T(n,n);
    
    C = A.triangularView<Lower>().solve(b);
    T = A.triangularView<StrictlyUpper>();
    T = A.triangularView<Lower>().solve(T);
    
    iter = 0 ;
    res = 1 ;
    while(res > tol && iter < max_iter){
           x = -T*(x) + C ;
           res = (x1-x).norm();
           x1 = x;
           iter++;

           if( iter % report_interval == 0 & report ) {
                       std::cout << "Iteration: " << iter << "     Residual: " << res << std::endl ;
           }
    }

    if( report ) {
      std::cout << "Iteration: " << iter << "     Residual: " << res << std::endl ;
    }

    if(res > tol) {
      std::cout << "=========================================================="<<std::endl ;
      std::cout << " Warning Solution not converged before maximum iterations "<<std::endl ;
      std::cout << "               Residual = "<< res <<std::endl ;
      std::cout << "=========================================================="<<std::endl ;
    }
}

void Solver::jacobi(){

    VectorXd x1 = VectorXd::Zero(n) ;
    VectorXd C(n) ;
    MatrixXd T(n,n);
    MatrixXd D(n,n);
    
    D = A.diagonal().asDiagonal();
    C = D.inverse() * b;
    T = D.inverse() * (A-D);
    
    iter = 0 ;
    res = 1 ;
    while(res > tol && iter < max_iter){
           x = -T*(x) + C ;
           res = (x1-x).norm();
           x1 = x;
           iter++;

           if( iter % report_interval == 0 & report ) {
                       std::cout << "Iteration: " << iter << "     Residual: " << res << std::endl ;
           } 

         }
    if( report ) {
      std::cout << "Iteration: " << iter << "     Residual: " << res << std::endl ;
    }

    if(res > tol) {
      std::cout << "=========================================================="<<std::endl ;
      std::cout << " Warning Solution not converged before maximum iterations "<<std::endl ;
      std::cout << "               Residual = "<< res <<std::endl ;
      std::cout << "=========================================================="<<std::endl ;
    }
}

void Solver::solution_init(int n_eq) {
    n = n_eq ;
    A = MatrixXd::Zero(n_eq, n_eq) ;
    b = VectorXd::Zero(n_eq) ;
    x = b ;
}

void Solver::solve() {
  switch(method){
    case 0:
      jacobi() ;
      break;
  case 1:
      gaussSeidel() ;
  break;
  }
}
