#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "elements.hpp"
#include "./eigen3/Eigen/Dense"
#include <vector>

TEST_CASE( "Node Default Constructor", "[Node]" ) {
    Node n(2) ;
    CHECK( n.ind == 2 ) ;
    CHECK( n.df == 0 ) ;
    CHECK( n.boundary == false ) ;
    CHECK( n.BC == 0 ) ;
}

TEST_CASE( "Node Full Constructor", "[Node]" ) {
    Node n1(5, 1, true, 4) ;
    CHECK( n1.ind == 5) ;
    CHECK( n1.df == 1 ) ;
    CHECK( n1.boundary == true ) ;
    CHECK( n1.BC == 4 ) ;
 }

TEST_CASE( "Node1d Default Constructor", "[Node][Node1d]" ) {
    Node1d n1d(0, 25.6) ;
    CHECK( n1d.ind == 0 ) ;
    CHECK( n1d.coords == Approx(25.6) ) ;
    CHECK( n1d.df == 0 ) ;
    CHECK( n1d.boundary == false ) ;
    CHECK( n1d.BC == 0 ) ;
}
TEST_CASE( "Node1d Full Constructor", "[Node][Node1d]" ) {
    Node1d n1d1(5, 0.000001, 1, true, 4) ;
    CHECK( n1d1.ind == 5 ) ;
    CHECK( n1d1.coords == Approx(0.000001) ) ;
    CHECK( n1d1.df == 1 ) ;
    CHECK( n1d1.boundary == true ) ;
    CHECK( n1d1.BC == 4 ) ;
}

TEST_CASE( "Node2d Default Constructor", "[Node][Node2d]" ) {
    Node2d n2d(0, 25.6, 17.1) ;
    CHECK( n2d.ind == 0 ) ;
    CHECK( n2d.coords(0) == Approx(25.6) ) ;
    CHECK( n2d.coords(1) == Approx(17.1) ) ;
    CHECK( n2d.df == 0 ) ;
    CHECK( n2d.boundary == false ) ;
    CHECK( n2d.BC == 0 ) ;
}

TEST_CASE( "Node2d Full Constructor", "[Node][Node2d]" ) {
    Node2d n2d1(5, 0.000001, 0.002, 1, true, 4) ;
    CHECK( n2d1.ind == 5 ) ;
    CHECK( n2d1.coords(0) == Approx(0.000001) ) ;
    CHECK( n2d1.df == 1 ) ;
    CHECK( n2d1.boundary == true ) ;
    CHECK( n2d1.BC == 4 ) ;
}

TEST_CASE( "Line Constructor with Node element", "[Line]") {
    Node1d n1(0, 0.0) ;
    Node1d n2(1, 0.1) ;
    Line1d L1(&n1, &n2) ;

    CHECK( L1.nodes[0] == &n1 ) ;
    CHECK( L1.nodes[0]->ind == n1.ind ) ;
    CHECK( L1.nodes[1] == &n2 ) ;
    CHECK( L1.nodes[1]->df == n2.df ) ;
}

TEST_CASE( "Test addNodes method of Line", "[Line]") {
    Node1d n1(0, 0.0) ;
    Node1d n2(1, 0.1) ;
    Node1d n3(2, 0.2) ;
    
    std::vector<Node1d> nodes_g ;
    nodes_g.push_back(n1);
    nodes_g.push_back(n2);
    nodes_g.push_back(n3);
    
    std::vector<Node1d*> nodes_e ;
    for(int i = 0 ; i<3 ; i++){
        nodes_e.push_back(&nodes_g[i]);
    }
    
    Line1d L1(&n1, &n2) ;
    
    L1.addNodes(2, nodes_g, nodes_e) ; 

    CHECK( L1.nodes.size() == 4 ) ;
    CHECK( nodes_e.size() == 5 ) ;
    CHECK( nodes_g.size() == 5 ) ;
}
        
        
TEST_CASE( "Test Instantiate an Element1d", "[Element1d]") {
    Node1d n0(0, 0.0) ;
    Node1d n1(1, 1.0) ;
    Node1d n2(2, 2.0) ;
    
    std::vector<Node1d*> nodes_g ;
    nodes_g.push_back(&n0);
    nodes_g.push_back(&n1);
    nodes_g.push_back(&n2);
    
    Line1d L0(&n0, &n1) ;
    Line1d L1(&n1, &n2) ;
    Line1d L2(&n2, &n0) ;

    std::vector<Line1d*> lines_g ;
    lines_g.push_back(&L0);
    lines_g.push_back(&L1);
    lines_g.push_back(&L2);
 
 	std::vector<Node1d*> nodes_e(&nodes_g[0],&nodes_g[2]) ;
 	std::vector<Line1d*> lines_e(&lines_g[0],&lines_g[1]) ;
    Element1d E(0, 1, 3, nodes_e, lines_e) ;  

    nodes_e.assign(&nodes_g[1],&nodes_g[3]) ;
 	lines_e.assign(&lines_g[1],&lines_g[2]) ;
    Element1d E1(1, 1, 3, nodes_e, lines_e) ; 

    CHECK( E.jac_det == 0.5 ) ;
    CHECK( E.Edges.size() == 1 ) ;
    CHECK( E.Nodes.size() == 2 ) ;

    CHECK( E.order == 1) ;
    CHECK( E.quad_pts == 3) ;
    
    E.quadrature(0) ;
    CHECK( E.w == Approx(0.8888888888888888888888889) );
    CHECK( E.xi == Approx(0.0) ) ;

    E.quadrature(1) ; 
    CHECK( E.xi == Approx(0.7745966692414833770358531) ) ;
    CHECK( E.w == Approx(0.5555555555555555555555556) ) ;

    MatrixXd A = MatrixXd::Zero(3,3) ;
    VectorXd b = VectorXd::Zero(3) ;

    E.AbCalc(A, b);
    MatrixXd A_e(3,3);
    VectorXd b_e(3);
    A_e << 1,-1,0,-1,1,0,0,0,0;
    b_e << 0.5,0.5,0;

    CHECK( (A-A_e).sum() == Approx(0.0) );
    CHECK( (b-b_e).sum() == Approx(0.0) );
    E1.AbCalc(A, b);

}

TEST_CASE( "Test Instantiate a Domain1d", "[Domain1d]") {
	
	int order_init = 1 ; 
	int nx_init = 3 ; 
	int quad_pts_init = 3 ; 
	double Xmin_init = 0.0 ; 
	double Xmax_init = 3.0 ; 
	double dirichlet_init = 0.0 ;

	Domain1d Domain(order_init, nx_init, quad_pts_init, Xmin_init, Xmax_init, dirichlet_init ) ;

	CHECK( Domain.Elements.capacity() == 3 ) ;
    CHECK( Domain.Nodes.capacity() == 4 ) ;
    CHECK( Domain.Edges.capacity() == 3 ) ;
	CHECK( Domain.Elements.size() == 0 ) ;
    CHECK( Domain.Nodes.size() == 0 ) ;
    CHECK( Domain.Edges.size() == 0 ) ;

    Domain.build_elements( ) ;
    CHECK( Domain.Elements.size() == 3 ) ;
    CHECK( Domain.Nodes.size() == 4 ) ;
    CHECK( Domain.Edges.size() == 3 ) ;

    MatrixXd A = MatrixXd::Zero(Domain.Nodes.size(), Domain.Nodes.size()) ;
    VectorXd b = VectorXd::Zero(Domain.Nodes.size()) ;

    Domain.build_Ab(&A, &b) ;

    MatrixXd A0(4,4) ;
    		 A0 << 1, -1,  0,  0,
	              -1,  2, -1,  0,
	               0, -1,  2, -1,
 	               0,  0, -1,  1;

	VectorXd b0(4) ;
			 b0 << 0.5, 1, 1, 0.5;

	CHECK( (A-A0).norm() == Approx(0.0) ) ;
    CHECK( (b-b0).norm() == Approx(0.0) ) ;

    Domain.add_constraints(A, b) ;

    MatrixXd A1(4,4) ;
    		 A1 << 1,  0,  0,  0,
	               0,  2, -1,  0,
	               0, -1,  2,  0,
 	               0,  0,  0,  1;

 	VectorXd b1(4) ;
 			 b1 << 0, 1, 1, 0 ;

 	CHECK( (A-A1).norm() == Approx(0.0) ) ;
 	CHECK( (b-b1).norm() == Approx(0.0) ) ;
    
}