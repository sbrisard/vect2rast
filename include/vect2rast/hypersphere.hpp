/**
 * @file v2r_ndsphere.h
 *
 * @brief This header defines n-dimensional spheres.
 */
#pragma once

#include <array>
#include <ostream>
#include <sstream>

#include "vect2rast.hpp"

template <size_t DIM>
class Hypersphere {
 public:
  /** Coordinates of the center. */
  const std::array<double, DIM> center;

  /** The radius of the hypersphere. */
  const double radius;

  Hypersphere(const std::array<double, DIM> center, const double radius)
      : center(center), radius(radius) {}

  std::string repr() const {
    std::ostringstream stream;
    stream << "Hypersphere<" << DIM << ">{radius=" << radius
           << ", center=" << vect2rast::repr(center.cbegin(), center.cend())
           << "}";
    return stream.str();
  }

  void get_bounding_box(double *bbmin, double *bbmax) const {
    for (size_t i = 0; i < DIM; i++) {
      bbmin[i] = center[i] - radius;
      bbmax[i] = center[i] + radius;
    }
  }

  bool belongs(double const *point) const {
    double r2 = 0.0;
    for (size_t i = 0; i < DIM; i++) {
      const double x_i = point[i] - center[i];
      r2 += x_i * x_i;
    }
    return r2 <= radius * radius;
  }
};

template <size_t DIM>
std::ostream &operator<<(std::ostream &os, const Hypersphere<DIM> hypersphere) {
  return os << hypersphere.repr();
}
