#define _USE_MATH_DEFINES

#include "catch2/catch.hpp"
#include "vect2rast/vect2rast.hpp"

TEST_CASE("Test case #1") {
  SECTION("Section #1") {
    REQUIRE(vect2rast::author() == "S. Brisard");
    REQUIRE(vect2rast::version() == "0.1");
    REQUIRE(vect2rast::return_one() == 1);
  }
}
