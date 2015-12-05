#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "elements.hpp"
#include "./eigen3/Eigen/Dense"
#include <vector>

TEST_CASE( "Node Instantiation", "[Node]" ) {
    Node n(2) ;
    CHECK( n.n == 2 ) ;
    CHECK( n.df == 0 ) ;
    CHECK( n.boundary == false ) ;
    CHECK( n.BC == 0 ) ;

    Node n1(5, 1, true, 4) ;
    CHECK( n1.n == 5) ;
    CHECK( n1.df == 1 ) ;
    CHECK( n1.boundary == true ) ;
    CHECK( n1.BC == 4 ) ;

 }

TEST_CASE( "Node1d Instantiation", "[Node]" ) {
    Node1d n1d(0, 25.6) ;
    CHECK( n1d.n == 0 ) ;
    CHECK( n1d.coords == Approx(25.6) ) ;
    CHECK( n1d.df == 0 ) ;
    CHECK( n1d.boundary == false ) ;
    CHECK( n1d.BC == 0 ) ;

    Node1d n1d1(5, 0.000001, 1, true, 4) ;
    CHECK( n1d1.n == 5 ) ;
    CHECK( n1d1.coords == Approx(0.000001) ) ;
    CHECK( n1d1.df == 1 ) ;
    CHECK( n1d1.boundary == true ) ;
    CHECK( n1d1.BC == 4 ) ;

 }
