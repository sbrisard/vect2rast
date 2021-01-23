#pragma once

#include <array>
#include <cmath>
#include <vector>

// DllExport void v2r_test_raster(V2R_Object const *object, double const
// *length,
//                               size_t const *size);

template <size_t DIM>
std::vector<std::array<double, DIM>> generate_directions() {
  if constexpr (DIM == 2) {
    size_t const num_directions = 10;
    std::vector<std::array<double, DIM>> directions{};
    for (size_t i = 0; i < num_directions; i++) {
      const double theta = 2 * M_PI * i / (double)num_directions;
      directions.push_back({cos(theta), sin(theta)});
    }
    return directions;
  } else if constexpr (DIM == 3) {
    const double phi = .5 * (1. + sqrt(5.));
    const double u = 1. / sqrt(1 + phi * phi);
    const double v = phi * u;
    return {{0., -u, -v}, {0., -u, +v}, {0., +u, -v}, {0., +u, +v},
            {-u, -v, 0.}, {-u, +v, 0.}, {+u, -v, 0.}, {+u, +v, 0.},
            {-v, 0., -u}, {-v, 0., +u}, {+v, 0., -u}, {+v, 0., +u}};
  } else {
    // This should never occur, since DIM was already statically asserted
    // through get_num_directions
    return nullptr;
  }
}

std::array<double, 3> v2r_cross(const std::array<double, 3>& v1,
                                const std::array<double, 3>& v2) {
  return {v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2],
          v1[0] * v2[1] - v1[1] * v2[0]};
}

void v2r_normalize(std::array<double, 3>& v) {
  double s = 1. / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

void assert_true(bool predicate) {
  if (!predicate) exit(-1);
}

void assert_false(bool predicate) {
  if (predicate) exit(-1);
}

void assert_equals_double(double expected, double actual, double rtol,
                          double atol) {
  if (fabs(actual - expected) > rtol * fabs(expected) + atol) {
    exit(-1);
  }
}