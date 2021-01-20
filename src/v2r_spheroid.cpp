#include <math.h>
#include <string.h>

#include "vect2rast/v2r_spheroid.hpp"

void *v2r_spheroid_data_new(double equatorial_radius, double polar_radius,
                            double const *axis) {
  size_t const dim = 3;
  size_t const size = dim * sizeof(double);
  auto data = static_cast<V2R_SpheroidData *>(malloc(sizeof(V2R_SpheroidData)));
  data->equatorial_radius = equatorial_radius;
  data->polar_radius = polar_radius;
  const double a2 = equatorial_radius * equatorial_radius;
  const double c2 = polar_radius * polar_radius;
  data->q1 = 1. / a2;
  data->q2 = 1. / c2 - data->q1;
  // TODO Normalize axis?
  data->axis = static_cast<double *>(malloc(size));
  memcpy(data->axis, axis, size);
  return data;
}

void *v2r_spheroid_data_copy(void const *data) {
  auto data_ = static_cast<V2R_SpheroidData const *>(data);
  auto copy = static_cast<V2R_SpheroidData *>(v2r_spheroid_data_new(
      data_->equatorial_radius, data_->polar_radius, data_->axis));
  return copy;
}

void v2r_spheroid_data_free(void *data) {
  auto spheroid_data = static_cast<V2R_SpheroidData *>(data);
  free(spheroid_data->axis);
  free(data);
}

void v2r_spheroid_get_bounding_box(V2R_Object const *spheroid, double *bbmin,
                                   double *bbmax) {
  const size_t dim = spheroid->type->dim;
  auto data = static_cast<V2R_SpheroidData const *>(spheroid->data);
  const double a = data->equatorial_radius;
  const double c = data->polar_radius;
  const double a2 = a * a;
  const double c2 = c * c;
  const double c2_m_a2 = c2 - a2;
  double *end = spheroid->center + dim;
  for (double *x = spheroid->center, *n = data->axis, *x1 = bbmin, *x2 = bbmax;
       x < end; x += 1, n += 1, x1 += 1, x2 += 1) {
    const double r = sqrt(a2 + c2_m_a2 * (*n) * (*n));
    *x1 = *x - r;
    *x2 = *x + r;
  }
}

static bool v2r_spheroid_belongs(V2R_Object const *spheroid,
                                 double const *point) {
  auto data = static_cast<V2R_SpheroidData const *>(spheroid->data);
  double x_dot_x = 0.;
  double d_dot_x = 0.;
  for (size_t i = 0; i < spheroid->type->dim; i++) {
    const double x_i = point[i] - spheroid->center[i];
    x_dot_x += x_i * x_i;
    d_dot_x += data->axis[i] * x_i;
  }
  return data->q1 * x_dot_x + data->q2 * d_dot_x * d_dot_x <= 1.;
}

double v2r_spheroid_equatorial_radius(V2R_Object const *spheroid) {
  return V2R_SPHEROID_DATA(spheroid)->equatorial_radius;
}

double v2r_spheroid_polar_radius(V2R_Object const *spheroid) {
  return V2R_SPHEROID_DATA(spheroid)->polar_radius;
}

void v2r_spheroid_axis(V2R_Object const *spheroid, double *axis) {
  V2R_SpheroidData *data = V2R_SPHEROID_DATA(spheroid);
  for (size_t i = 0; i < spheroid->type->dim; i++) {
    axis[i] = data->axis[i];
  }
}

V2R_ObjectType Spheroid = {.name = "Spheroid",
                           .dim = 3,
                           .data_copy = v2r_spheroid_data_copy,
                           .data_free = v2r_spheroid_data_free,
                           .belongs = v2r_spheroid_belongs,
                           .get_bounding_box = v2r_spheroid_get_bounding_box};

V2R_Object *v2r_spheroid_new(double const *center, double equatorial_radius,
                             double polar_radius, double const *axis) {
  V2R_Object *object = v2r_object_new(&Spheroid);
  object->data = v2r_spheroid_data_new(equatorial_radius, polar_radius, axis);
  memcpy(object->center, center, Spheroid.dim * sizeof(double));
  return object;
}
