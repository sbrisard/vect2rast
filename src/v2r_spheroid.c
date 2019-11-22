#include <math.h>

#include "vect2rast.h"

#define V2R_SPHEROID_DATA(spheroid) ((V2R_SpheroidData *)((spheroid)->data))

typedef struct V2R_SpheroidData_ {
  double equatorial_radius;
  double polar_radius;
  double *axis;
  double q1, q2;
} V2R_SpheroidData;

void v2r_spheroid_data_free(void *data) {
  V2R_SpheroidData *spheroid_data = data;
  free(spheroid_data->axis);
  free(data);
}

static bool v2r_spheroid_belongs(V2R_Object const *spheroid,
                                 double const *point) {
  V2R_SpheroidData const *data = spheroid->data;
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

V2R_ObjectType const Spheroid = {.name = "Spheroid",
                                 .dim = 3,
                                 .belongs = v2r_spheroid_belongs,
                                 .data_free = v2r_spheroid_data_free};

V2R_Object *v2r_spheroid_new(double const *center, double equatorial_radius,
                             double polar_radius, double const *axis) {
  V2R_SpheroidData *data = malloc(sizeof(V2R_SpheroidData));
  V2R_Object *object = v2r_object_new(&Spheroid);
  data->axis = malloc(Spheroid.dim * sizeof(double));
  object->data = data;

  data->equatorial_radius = equatorial_radius;
  data->polar_radius = polar_radius;
  const double a2 = equatorial_radius * equatorial_radius;
  const double c2 = polar_radius * polar_radius;
  const double c2_m_a2 = c2 - a2;
  data->q1 = 1. / a2;
  data->q2 = 1. / c2 - data->q1;
  for (size_t i = 0; i < Spheroid.dim; i++) {
    /* TODO Normalize axis? */
    data->axis[i] = axis[i];
    object->center[i] = center[i];
    const double r = sqrt(a2 + c2_m_a2 * axis[i] * axis[i]);
    object->bbmin[i] = object->center[i] - r;
    object->bbmax[i] = object->center[i] + r;
  }
  return object;
}
