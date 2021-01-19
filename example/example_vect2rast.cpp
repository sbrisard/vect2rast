#include <iostream>
#include "vect2rast/vect2rast.hpp"

int main() {
  std::cout << "version = " << vect2rast::version() << std::endl;
  std::cout << "author = " << vect2rast::author() << std::endl;
  std::cout << "return_one() = " << vect2rast::return_one() << std::endl;
  return 0;
}
