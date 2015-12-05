#define CATCH_CONFIG_MAIN
#include "./external/catch.hpp"
#include "factorial.hpp"

TEST_CASE( "Factorials are computed", "[factorial]") {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}


TEST_CASE( "More Factorials", "[factorial]") {
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}
