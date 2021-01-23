/**
 * @file v2r_spheroid.h
 *
 * @brief This header defines n-dimensional spheres.
 */
#ifndef __V2R_SPHEROID_H_202001120654__
#define __V2R_SPHEROID_H_202001120654__

#include <array>
#include <span>

#include "vect2rast.hpp"

#define V2R_SPHEROID_DATA(spheroid) ((V2R_SpheroidData *)((spheroid)->data))

typedef struct V2R_SpheroidData_ {
  double equatorial_radius;
  double polar_radius;
  double *axis;
  double q1, q2;
} V2R_SpheroidData;

DllExport double v2r_spheroid_equatorial_radius(V2R_Object const *spheroid);
DllExport double v2r_spheroid_polar_radius(V2R_Object const *spheroid);
DllExport void v2r_spheroid_axis(V2R_Object const *spheroid, double *axis);
DllExport V2R_Object *v2r_spheroid_new(double const *center,
                                       double equatorial_radius,
                                       double polar_radius, double const *axis);
namespace vect2rast {
class Spheroid {
 public:
  const std::array<double, 3> center;
  const double equatorial_radius;
  const double polar_radius;
  const std::array<double, 3> axis;

  Spheroid(const std::array<double, 3> center, double equatorial_radius,
           double polar_radius, const std::array<double, 3> axis)
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
  void get_bounding_box(std::span<double, 3> bbmin,
                        std::span<double, 3> bbmax) const {
    const size_t dim = 3;
    const double a2 = equatorial_radius * equatorial_radius;
    const double c2 = polar_radius * polar_radius;
    const double c2_m_a2 = c2 - a2;
    for (size_t i = 0; i < dim; i++) {
      const double r = sqrt(a2 + c2_m_a2 * axis[i] * axis[i]);
      bbmin[i] = center[i] - r;
      bbmax[i] = center[i] + r;
    }
  }

  /** Return `true` if the specified `point` belongs to this hypersphere. */
  bool belongs(const std::span<double, 3> point) const {
    double x_dot_x = 0.;
    double d_dot_x = 0.;
    for (size_t i = 0; i < 3; i++) {
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
}  // namespace vect2rast
#endif
