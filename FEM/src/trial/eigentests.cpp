#include "./external/eigen3/Eigen/Dense"
#include <cmath>
#include <iostream>

using namespace Eigen;

int gaussSeidel(MatrixXd& A, VectorXd& b, VectorXd& x0);
int jacobi(MatrixXd& A, VectorXd& b, VectorXd& x0);

int main(){
    MatrixXd m1(4,4);
    int n = m1.cols() ;
    VectorXd b(n) ;
    VectorXd xj = VectorXd::Zero(n) ;
    VectorXd xg = VectorXd::Zero(n) ;
    double res = 1 ;
    double limit = 1e-11;
    
    m1 << 10, -1, 2, 0, -1 , 11, -1, 3, 2, -1, 10, -1, 0, 3, -1, 8;
    b << 6, 25, -11, 15;
    
    std::cout << "===================Gauss-Seidel=====================" << std::endl;
    gaussSeidel(m1, b, xg) ;    
    std::cout << "======================Jacobi========================" << std::endl;
    jacobi(m1, b, xj) ;  
}
    
int gaussSeidel(MatrixXd &A, VectorXd& b, VectorXd& x0){
    int n = A.cols() ;
    VectorXd x1 = VectorXd::Zero(n) ;
    VectorXd C(n) ;
    MatrixXd T(n,n);

    double res = 1 ;
    double tol = 1e-11;
    int max_iter = 20 ;
    
    C = A.triangularView<Lower>().solve(b);
    T = A.triangularView<StrictlyUpper>();
    T = A.triangularView<Lower>().solve(T);
    
    std::cout << "T = " << T << std::endl;
    std::cout << "C = " << C << std::endl;
    
    int iter = 0 ;
    while(res > tol && iter < max_iter){
           x0 = -T*(x0) + C ;
           res = (x1-x0).array().abs().sum();
           x1 = x0;
           
           std::cout << "residual = " << res << std::endl;
           std::cout << "x = " << x1 << std::endl;
           iter++;
    }
}

int jacobi(MatrixXd &A, VectorXd& b, VectorXd& x0){
    int n = A.cols() ;
    VectorXd x1 = VectorXd::Zero(n) ;
    VectorXd C(n) ;
    MatrixXd T(n,n);
    MatrixXd D(n,n);
    double res = 1 ;
    double tol = 1e-11;
    int max_iter = 20 ;
    
    D = A.diagonal().asDiagonal();
    C = D.inverse()*b;
    T = D.inverse()*(A-D);
 
    std::cout << "T = \n" << T << std::endl;
    std::cout << "C = \n" << C << std::endl;
    
    int iter = 0 ;
    while(res > tol && iter < max_iter){
           x0 = -T*(x0) + C ;
           res = (x1-x0).array().abs().sum();
           x1 = x0;
           
           std::cout << "residual = " << res << std::endl;
           std::cout << "x = \n" << x1 << std::endl;
           iter++;
    }
}
