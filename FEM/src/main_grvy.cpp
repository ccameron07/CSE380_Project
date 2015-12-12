#include "elements.hpp"
#include "solvers.hpp"
#include "masa_helper.hpp"
#include <./eigen3/Eigen/Dense>
#include <masa.h>
#include <grvy.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define FUNC_BEGIN_TIMER gt.BeginTimer(__func__);
#define FUNC_END_TIMER   gt.EndTimer  (__func__);

using namespace MASA ;
using namespace GRVY ;

GRVY::GRVY_Timer_Class gt;
void setup (Domain1d &D, Solver &S) ;

int main() {

    gt.Init("1d FEM Timer");

    gt.BeginTimer("Main Program");

    gt.BeginTimer("Initialization");
    Domain1d Domain ;
    Solver Solution ;
    
    setup(Domain, Solution) ;
   
    Domain.Elements.reserve(Domain.nx) ;
    Domain.Nodes.reserve(Domain.order*Domain.nx+1) ;
    Domain.Edges.reserve(Domain.nx) ;

    if(Domain.with_masa){
        Manufactured mfg("Chris Test 2","heateq_1d_steady_const");
        Domain.stiffness = [&mfg] (double x)->double{return mfg.stiffness;} ;
        Domain.forcing = [&mfg] (double x)->double{return mfg.forcing(x);} ;
    } else {
        Domain.stiffness = [] (double x)->double{return 1;} ;
        Domain.forcing = [] (double x)->double{return 1;} ;
    }
    gt.EndTimer("Initialization") ;
    
    gt.BeginTimer("Build Elements") ;
    Domain.build_elements( ) ;
    if (Domain.order == 2) {
        Domain.addNodes(1);
    }
    gt.EndTimer("Build Elements") ;
    
    gt.BeginTimer("Assemble Equations");
 	Solution.solution_init(Domain.Nodes.size()) ;
	Domain.build_Ab(&Solution.A, &Solution.b);
	Domain.add_constraints(Solution.A, Solution.b);
    std::cout << Solution.A << std::endl ;
    gt.EndTimer("Assemble Equations");

    gt.BeginTimer("Solve Equations");
 	Solution.solve() ;
    gt.EndTimer("Solve Equations");

    Solution.output(Domain);
    gt.EndTimer("Main Program");
    gt.Finalize();
    gt.Summarize();
    gt.Reset();

    return 0;
}


void setup (Domain1d &D, Solver &S){
	
    GRVY_Input_Class iparse;
   
    if(! iparse.Open("./input_file.txt"))
          exit(1);

    if ( iparse.Read_Var("with_masa",&D.with_masa,true) )
        std::cout << "Use MASA?  = " << D.with_masa << std::endl ;

    S.with_masa = D.with_masa ;

    if ( iparse.Read_Var("output_file",&S.file_out,"output.txt") )
        std::cout << "Output File = " << S.file_out << std::endl ;
    
    if ( iparse.Read_Var("n_x",&D.nx,8) )
        std::cout << "n_x = " << D.nx << std::endl ;
    
    if ( iparse.Read_Var("Xmin",&D.Xmin) )
        std::cout << "Xmin = " << D.Xmin << std::endl ;

    if ( iparse.Read_Var("Xmax",&D.Xmax) )
        std::cout << "Xmax = " << D.Xmax << std::endl ;
    
    if ( iparse.Read_Var("Dirichlet",&D.dirichlet,0.0) )
        std::cout << "Dirichlet BC  = " << D.dirichlet << std::endl ;

    if ( iparse.Read_Var("order",&D.order,2) )
        std::cout << "Shape function order =  " << D.order << std::endl ;

    if ( iparse.Read_Var("quad",&D.quad_pts,8) )
        std::cout << "Number of quadrature points = " << D.quad_pts << std::endl ;
    
    if ( iparse.Read_Var("method",&S.method,8) )
        std::cout << "Solution method = " << S.method << std::endl ;
    
    if ( iparse.Read_Var("tol",&S.tol,1e-14) )
        std::cout << "Tolerance (methods 1,2) = " << S.tol << std::endl ;
    
    if ( iparse.Read_Var("max_iter",&S.max_iter) )
        std::cout << "Solution maximum iterations = " << S.max_iter << std::endl ;
    
    if ( iparse.Read_Var("report",&S.report,true) )
        std::cout << "Report convergence = " << S.report << std::endl ;
    
    if ( iparse.Read_Var("report_interval",&S.report_interval,1000) )
        std::cout << "Report Interval = " << S.report_interval << std::endl ;
}
    /*
    int order_init = 2 ; 
	int nx_init = 8; 
	int quad_pts_init = 4 ; 
	double Xmin_init = 0 ; 
	double Xmax_init = 4 ; 
	double dirichlet_init = 0.0 ;
	int method_init = 3 ;
	double tol_init = 1e-14 ;
	int max_iter_init = 1000000 ;
    int report_interval = 1000 ;
    bool report = true ;
    bool with_masa = true ;
    std::string file_out = "masa_validate_2nd_sparse_8.txt";
    */












