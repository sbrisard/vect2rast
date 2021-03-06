#pragma once

#include <array>
#include <cmath>
#include <numeric>
#include <vector>

#include "vect2rast/vect2rast.hpp"

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

template <size_t N>
void v2r_normalize(std::array<double, N>& v) {
  double s =
      1. / sqrt(std::inner_product(v.cbegin(), v.cend(), v.cbegin(), 0.0));
  std::transform(v.cbegin(), v.cend(), v.begin(),
                 [s](double x) { return s * x; });
}
