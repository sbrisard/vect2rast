#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "vect2rast.h"
#include "v2r_test_utils.h"

void v2r_test_spheroid_new() {
  size_t const dim = 3;
  double const x[] = {1.2, -3.4, 5.6};
  double const theta = .3 * M_PI;
  double const phi = .5 * M_PI;
  double const n[] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
  double const a = 7.8;
  double const c = 9.1;

  V2R_Object *spheroid = v2r_spheroid_new(x, a, c, n);
  g_assert_cmpint(spheroid->type->dim, ==, dim);
  g_assert_cmpfloat(v2r_spheroid_equatorial_radius(spheroid), ==, a);
  g_assert_cmpfloat(v2r_spheroid_polar_radius(spheroid), ==, c);

  double n_act[3];
  v2r_spheroid_axis(spheroid, n_act);
  for (size_t i = 0; i < dim; i++) {
    g_assert_cmpfloat(spheroid->center[i], ==, x[i]);
    g_assert_cmpfloat(n_act[i], ==, n[i]);
  }

  v2r_object_free(spheroid);
}

void v2r_setup_test_spheroid_belongs() {
  double center[] = {1.2, -3.4, 5.6};
  double axis[] = {1., 0., 0.};
  V2R_Object *spheroid = v2r_spheroid_new(center, 1.0, 0.1, axis);
  void *data = v2r_test_belongs_data_new(spheroid, center, true);
  g_test_add_data_func_full("/spheroid/belongs/1", data, v2r_test_belongs,
                            v2r_test_belongs_data_free);
}

void v2r_setup_test_spheroid() {
  g_test_add_func("/spheroid/new", v2r_test_spheroid_new);
  v2r_setup_test_spheroid_belongs();
}
