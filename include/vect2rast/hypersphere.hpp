/** N-dimensional sphere. */
#pragma once

#include <array>
#include <numeric>
#include <ostream>
#include <span>
#include <sstream>

#include "vect2rast.hpp"

namespace vect2rast {
template <size_t DIM>
class Hypersphere {
 public:
  /** Number of spatial dimensions. */
  static constexpr size_t dim = DIM;

  /** Coordinates of the center. */
  const std::array<double, DIM> center;

  /** The radius of the hypersphere. */
  const double radius;

  /**
   * Create a new instance of this class, with specified `center` and `radius`.
   */
  Hypersphere(const std::array<double, DIM> center, const double radius)
      : center(center), radius(radius) {}

  /** Return a string representation of this hypersphere. */
  std::string repr() const {
    std::ostringstream stream;
    stream << "Hypersphere<" << DIM << ">{radius=" << radius
           << ", center=" << vect2rast::repr(center.cbegin(), center.cend())
           << "}";
    return stream.str();
  }

  /** Updates `bbmin` and `bbmax` with the bounding-box of this hypersphere. */
  void get_bounding_box(std::span<double, DIM> bbmin,
                        std::span<double, DIM> bbmax) const {
    std::transform(center.cbegin(), center.cend(), bbmin.begin(),
                   [this](double x) { return x - this->radius; });
    std::transform(center.cbegin(), center.cend(), bbmax.begin(),
                   [this](double x) { return x + this->radius; });
  }

  /** Return `true` if the specified `point` belongs to this hypersphere. */
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
}  // namespace vect2rast
