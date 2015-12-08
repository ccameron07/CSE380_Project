#include "elements.hpp"
#include "solvers.hpp"
#include "masa_helper.hpp"
#include <./eigen3/Eigen/Dense>
#include <masa.h>
#include <grvy.h>
#include <iostream>
#include <vector>
#include<sys/time.h>
#include<time.h>
#include<unistd.h>

using namespace MASA ;
using namespace GRVY ;

#define FUNC_BEGIN_TIMER gt.BeginTimer(__func__);
#define FUNC_END_TIMER   gt.EndTimer  (__func__);

int main() {
    GRVY::GRVY_Timer_Class gt;

    //grvy_log_setlevel(GRVY_INFO);
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
    
    gt.Init("FEM 1D");    

    gt.BeginTimer("Main Program");
    gt.BeginTimer("Initialization");

    grvy_printf(GRVY_INFO, "Heyo does this break everything??\n") ;
    Manufactured mfg("MASA Solution","heateq_1d_steady_const");
    std::function<double(double)> stiffness = [&mfg] (double x)->double{return mfg.stiffness;} ;
    std::function<double(double)> forcing = [&mfg] (double x)->double{return mfg.forcing(x);} ;
	
    grvy_printf(GRVY_INFO, "Heyo does this break everything??\n") ;
    Domain1d Domain(order_init, nx_init, quad_pts_init, Xmin_init, Xmax_init, dirichlet_init ) ;
	Domain.forcing = forcing ;
	Domain.stiffness = stiffness ; 
    
    grvy_printf(GRVY_INFO, "Heyo does this break everything??\n") ;
    Domain.build_elements( ) ;

	Solver GaussSeidel(method_init, tol_init, max_iter_init, report_interval, report) ;
 	GaussSeidel.solution_init(Domain.Nodes.size()) ;

    grvy_printf(GRVY_INFO, "Heyo does this break everything??\n") ;
	Domain.build_Ab(&GaussSeidel.A, &GaussSeidel.b);
	Domain.add_constraints(GaussSeidel.A, GaussSeidel.b);
    gt.EndTimer("Initialization");

    grvy_printf(GRVY_INFO, "Heyo does this break everything??\n") ;
    gt.BeginTimer("Solver");
 	GaussSeidel.solve() ;
    gt.EndTimer("Solver");
 
    std::cout << GaussSeidel.x << std::endl;
    std::cout << "-------------------" << std::endl;
    mfg.exact_solution(Domain.Nodes) ;
    for(int i = 0; i < Domain.Nodes.size(); i++){
        std::cout << Domain.Nodes[i].coords << std::endl ;
    
    gt.EndTimer("Main Program");
    gt.Finalize();
    gt.Summarize();

    }
}
