#pragma once

#include <array>
#include <cmath>
#include <vector>
#include "vect2rast/vect2rast.hpp"

// DllExport void v2r_test_raster(V2R_Object const *object, double const
// *length,
//                               size_t const *size);

template <size_t DIM>
constexpr size_t get_num_directions() {
  static_assert((DIM == 2) || (DIM == 3));
  if constexpr (DIM == 2) {
    return 10;
  } else if constexpr (DIM == 3) {
    return 12;
  } else {
    // This should never occur
  }
}

template <size_t DIM>
std::vector<std::array<double, DIM>> generate_directions() {
  size_t const num_directions = get_num_directions<DIM>();
  std::vector<std::array<double, DIM>> directions{};
  if constexpr (DIM == 2) {
    for (size_t i = 0; i < num_directions; i++) {
      const double theta = 2 * M_PI * i / (double)num_directions;
      directions.push_back({cos(theta), sin(theta)});
    }
  } else if constexpr (DIM == 3) {
    const double phi = .5 * (1. + sqrt(5.));
    const double u = 1. / sqrt(1 + phi * phi);
    const double v = phi * u;
    directions.push_back({0., -u, -v});
    directions.push_back({0., -u, +v});
    directions.push_back({0., +u, -v});
    directions.push_back({0., +u, +v});
    directions.push_back({-u, -v, 0.});
    directions.push_back({-u, +v, 0.});
    directions.push_back({+u, -v, 0.});
    directions.push_back({+u, +v, 0.});
    directions.push_back({-v, 0., -u});
    directions.push_back({-v, 0., +u});
    directions.push_back({+v, 0., -u});
    directions.push_back({+v, 0., +u});
  } else {
    // This should never occur, since DIM was already statically asserted
    // through get_num_directions
    return nullptr;
  }
  return directions;
}

DllExport void v2r_cross(const std::array<double, 3>& v1,
                         const std::array<double, 3>& v2,
                         std::array<double, 3>& v3);
DllExport void v2r_normalize(std::array<double, 3>& v);

DllExport void assert_true(bool predicate);
DllExport void assert_false(bool predicate);
DllExport void assert_equals_double(double expected, double actual, double rtol,
                                    double atol);
