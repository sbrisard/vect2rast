#pragma once

#include <stdio.h>
#include <string.h>
#include <cmath>
#include <span>
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
double *generate_directions() {
  size_t const num_directions = get_num_directions<DIM>();
  size_t const size = DIM * num_directions * sizeof(double);
  auto directions = static_cast<double *>(malloc(size));
  if constexpr (DIM == 2) {
    double *dir = directions;
    for (size_t i = 0; i < num_directions; i++) {
      double const theta = 2 * M_PI * i / (double)num_directions;
      *dir = cos(theta);
      dir += 1;
      *dir = sin(theta);
      dir += 1;
    }
  } else if constexpr (DIM == 3) {
    double phi = .5 * (1. + sqrt(5.));
    double u = 1. / sqrt(1 + phi * phi);
    double v = phi * u;
    double dir[] = {0., -u, -v, 0., -u, +v, 0., +u, -v, 0., +u, +v,
                    -u, -v, 0., -u, +v, 0., +u, -v, 0., +u, +v, 0.,
                    -v, 0., -u, -v, 0., +u, +v, 0., -u, +v, 0., +u};
    memcpy(directions, dir, size);
  } else {
    // This should never occur, since DIM was already statically asserted
    // through get_num_directions
    return nullptr;
  }
  return directions;
}

DllExport void v2r_cross(const std::span<double, 3> v1,
                         const std::span<double, 3> v2,
                         std::span<double, 3> v3);
DllExport void v2r_normalize(std::span<double, 3> v);

DllExport void assert_true(bool predicate);
DllExport void assert_false(bool predicate);
DllExport void assert_equals_double(double expected, double actual, double rtol,
                                    double atol);
