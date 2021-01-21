/**
 * @file v2r_ndsphere.h
 *
 * @brief This header defines n-dimensional spheres.
 */
#pragma once

#include <array>
#include <numeric>
#include <ostream>
#include <span>
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

  void get_bounding_box(std::span<double, DIM> bbmin,
                        std::span<double, DIM> bbmax) const {
    std::transform(center.cbegin(), center.cend(), bbmin.begin(),
                   [=](double x) -> double { return x - radius; });
    std::transform(center.cbegin(), center.cend(), bbmax.begin(),
                   [=](double x) -> double { return x + radius; });
  }

  bool belongs(const std::span<double, DIM> point) const {
    return std::transform_reduce(point.begin(), point.end(), center.cbegin(),
                                 0.0, std::plus<>(), [](double c, double p) {
                                   double cp = p - c;
                                   return cp * cp;
                                 }) <= radius * radius;
  }
};

template <size_t DIM>
std::ostream &operator<<(std::ostream &os, const Hypersphere<DIM> hypersphere) {
  return os << hypersphere.repr();
}
