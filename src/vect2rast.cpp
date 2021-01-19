#include "vect2rast/vect2rast.hpp"

namespace vect2rast
{
  std::string author() {
    return PROJECT_AUTHOR_QUOTED;
  }

  std::string version() {
    // This preprocessor macro is defined in the CMakeLists.txt file
    return PROJECT_VERSION_QUOTED;
  }

  int return_one() { return 1; }
}  // namespace vect2rast
