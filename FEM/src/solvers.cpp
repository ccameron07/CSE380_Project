#include "elements.hpp"
#include "./eigen3/Eigen/Dense"
#include "./eigen3/Eigen/Sparse"
#include <iostream>
#include "solvers.hpp"
#include <grvy.h>
#include <masa.h>
#include <fstream>

using namespace MASA;
using namespace Eigen;
using namespace GRVY;

//Constuctor method
Solver::Solver(int method_init, double tol_init, int max_iter_init, int report_interval_init, bool report_init) {
  method = method_init ;
  tol = tol_init ;
  max_iter = max_iter_init ;
  report_interval = report_interval_init ;
  report = report_init ;
}

//Call Eigen Householder solution routine
void Solver::Householder(){

    x = A.colPivHouseholderQr().solve(b) ;
}

//Call Eigen conjugate gradient solver
void Solver::CG(){
    ConjugateGradient<SparseMatrix<double>> CG;

    CG.compute(A.sparseView());
    if(CG.info()!=Success) {
        std::cout << "decomposition failed" << std::endl ;
        return;
     }
     x = CG.solve(b);
     if(CG.info()!=Success) {
        std::cout << "solving failed" << std::endl ;
      return;
     }
}


void Solver::gaussSeidel(){

    VectorXd x1 = VectorXd::Zero(n) ;
    VectorXd C(n) ;
    MatrixXd T(n,n);
    SparseMatrix<double> Ts;

    //Create T matrix and C vector for performing terations
    C = A.triangularView<Lower>().solve(b);
    T = A.triangularView<StrictlyUpper>();
    T = A.triangularView<Lower>().solve(T);
    Ts = T.sparseView(); 
    
    iter = 0 ;
    res = 1 ;
    //iterate until tolerance max iterations reached
    while(res > tol && iter < max_iter){
           x = -Ts*(x) + C ;
           res = (x1-x).norm();
           x1 = x;
           iter++;

           //Intermittent convergence reporting
           if( iter % report_interval == 0 & report ) {
                       std::cout << "Iteration: " << iter << "     Residual: " << res << std::endl ;
           }
    }

    //Final convergence reporting
    if( report ) {
      std::cout << "Iteration: " << iter << "     Residual: " << res << std::endl ;
    }

    //Warning if tolerance not reached
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
    
    //Create Iterator matrix and vector
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

           //Intermittent Convergence reporting
           if( iter % report_interval == 0 & report ) {
                       std::cout << "Iteration: " << iter << "     Residual: " << res << std::endl ;
           } 

         }
    //Final convergence report     
    if( report ) {
      std::cout << "Iteration: " << iter << "     Residual: " << res << std::endl ;
    }

    //Warning if tolerance not reached
    if(res > tol) {
      std::cout << "=========================================================="<<std::endl ;
      std::cout << " Warning Solution not converged before maximum iterations "<<std::endl ;
      std::cout << "               Residual = "<< res <<std::endl ;
      std::cout << "=========================================================="<<std::endl ;
    }
}

//Method to initialize A, x, and b and reserve appropriate memory
void Solver::solution_init(int n_eq) {
    n = n_eq ;
    A = MatrixXd::Zero(n_eq, n_eq) ;
    b = VectorXd::Zero(n_eq) ;
    x = b ;
}

//Call the solver indicated by method
void Solver::solve() {
  switch(method){
    case 0:
      jacobi() ;
      break;
    case 1:
      gaussSeidel() ;
      break;
    case 2:
      Householder() ;
      break;
    case 3:
      Householder() ;
      break;
  }
}

//Reads in domain and outputs solution and exact solution (if using MASA)
//Also reports coordinates of each point in the domain
void Solver::output(Domain1d &Domain){
    
    std::ofstream f ;
    f.open(file_out);
    f.precision(16) ;
    
    if(with_masa){
        VectorXd exact = VectorXd::Zero(Domain.Nodes.size()) ;
        for (int i = 0; i < Domain.Nodes.size(); i++) {
            exact(i) = masa_eval_1d_exact_t( Domain.Nodes[i].coords ) ;
        }
    
        f << "x,exact,approx" << std::endl ;
        for(int i = 0; i < Domain.Nodes.size(); i++){
            f <<  Domain.Nodes[i].coords << "," << exact(i) << ',' << x(i) << std::endl ;
        }
    } else {
        f << "x,T" << std::endl ;
        for(int i = 0; i < Domain.Nodes.size(); i++){
        f <<  Domain.Nodes[i].coords << "," << x(i) << std::endl ;
        }
    }
    f.close();
}
