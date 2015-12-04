#include "catch.hpp"
#include "factorial.hpp"

TEST_CASE( "Factorial of zero is 0", "[factorial]") {
    REQUIRE( Factorial(0) == 1 );
}

