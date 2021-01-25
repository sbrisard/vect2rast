/** Define a 3D spheroid. */
#pragma once
#include <array>
#include <span>

#include "vect2rast.hpp"

namespace vect2rast {
class Spheroid {
 public:
  static constexpr size_t dim = 3;
  const std::array<double, Spheroid::dim> center;
  const double equatorial_radius;
  const double polar_radius;
  const std::array<double, Spheroid::dim> axis;

  Spheroid(const std::array<double, Spheroid::dim> center,
           double equatorial_radius, double polar_radius,
           const std::array<double, Spheroid::dim> axis)
      : center(center),
        equatorial_radius(equatorial_radius),
        polar_radius(polar_radius),
        axis(axis),
        q1(1. / (equatorial_radius * equatorial_radius)),
        q2(1. / (polar_radius * polar_radius) - q1) {}

  /** Return a string representation of this spheroid. */
  std::string repr() const {
    std::ostringstream stream;
    stream << "Spheroid"
           << "{center =" << vect2rast::repr(center.cbegin(), center.cend())
           << ", equatorial_radius =" << equatorial_radius
           << ", polar_radius =" << polar_radius
           << ", axis =" << vect2rast::repr(axis.cbegin(), axis.cend()) << "}";
    return stream.str();
  }

  /** Updates `bbmin` and `bbmax` with the bounding-box of this hypersphere. */
  void get_bounding_box(std::span<double, Spheroid::dim> bbmin,
                        std::span<double, Spheroid::dim> bbmax) const {
    const double a2 = equatorial_radius * equatorial_radius;
    const double c2 = polar_radius * polar_radius;
    const double c2_m_a2 = c2 - a2;
    for (size_t i = 0; i < Spheroid::dim; i++) {
      const double r = sqrt(a2 + c2_m_a2 * axis[i] * axis[i]);
      bbmin[i] = center[i] - r;
      bbmax[i] = center[i] + r;
    }
  }

  /** Return `true` if the specified `point` test_hypersphere_belongs to this hypersphere. */
  bool belongs(const std::span<double, Spheroid::dim> point) const {
    double x_dot_x = 0.;
    double d_dot_x = 0.;
    for (size_t i = 0; i < Spheroid::dim; i++) {
      const double x_i = point[i] - center[i];
      x_dot_x += x_i * x_i;
      d_dot_x += axis[i] * x_i;
    }
    return q1 * x_dot_x + q2 * d_dot_x * d_dot_x <= 1.;
  }

 private:
  const double q1;
  const double q2;
};

template <size_t DIM>
std::ostream &operator<<(std::ostream &os, const Spheroid spheroid) {
  return os << spheroid.repr();
}
}  // namespace vect2rast
