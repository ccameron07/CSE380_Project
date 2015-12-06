#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "elements.hpp"
#include "./eigen3/Eigen/Dense"
#include <vector>

TEST_CASE( "Node Default Constructor", "[Node]" ) {
    Node n(2) ;
    CHECK( n.n == 2 ) ;
    CHECK( n.df == 0 ) ;
    CHECK( n.boundary == false ) ;
    CHECK( n.BC == 0 ) ;
}

TEST_CASE( "Node Full Constructor", "[Node]" ) {
    Node n1(5, 1, true, 4) ;
    CHECK( n1.n == 5) ;
    CHECK( n1.df == 1 ) ;
    CHECK( n1.boundary == true ) ;
    CHECK( n1.BC == 4 ) ;
 }

TEST_CASE( "Node1d Default Constructor", "[Node][Node1d]" ) {
    Node1d n1d(0, 25.6) ;
    CHECK( n1d.n == 0 ) ;
    CHECK( n1d.coords == Approx(25.6) ) ;
    CHECK( n1d.df == 0 ) ;
    CHECK( n1d.boundary == false ) ;
    CHECK( n1d.BC == 0 ) ;
}
TEST_CASE( "Node1d Full Constructor", "[Node][Node1d]" ) {
    Node1d n1d1(5, 0.000001, 1, true, 4) ;
    CHECK( n1d1.n == 5 ) ;
    CHECK( n1d1.coords == Approx(0.000001) ) ;
    CHECK( n1d1.df == 1 ) ;
    CHECK( n1d1.boundary == true ) ;
    CHECK( n1d1.BC == 4 ) ;
}

TEST_CASE( "Node2d Default Constructor", "[Node][Node2d]" ) {
    Node2d n2d(0, 25.6, 17.1) ;
    CHECK( n2d.n == 0 ) ;
    CHECK( n2d.coords(0) == Approx(25.6) ) ;
    CHECK( n2d.coords(1) == Approx(17.1) ) ;
    CHECK( n2d.df == 0 ) ;
    CHECK( n2d.boundary == false ) ;
    CHECK( n2d.BC == 0 ) ;
}

TEST_CASE( "Node2d Full Constructor", "[Node][Node2d]" ) {
    Node2d n2d1(5, 0.000001, 0.002, 1, true, 4) ;
    CHECK( n2d1.n == 5 ) ;
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
    CHECK( L1.nodes[0]->n == n1.n ) ;
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
