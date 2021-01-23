#pragma once

#include <iterator>
#include <sstream>
#include <string>

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

// template <typename T>
// int v2r_raster(const T &object, double const *length, size_t const *size,
//               int *grid, int value);

namespace vect2rast {
template <typename Iterator>
std::string repr(Iterator first, Iterator last) {
  std::ostringstream stream;
  stream << "{";
  for (auto i = first; i < last; i++) {
    stream << *i << ",";
  }
  stream << "}";
  return stream.str();
}
}  // namespace vect2rast
