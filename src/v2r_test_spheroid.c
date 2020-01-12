#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "v2r_test_utils.h"
#include "v2r_spheroid.h"

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

void v2r_setup_test_spheroid_belongs(double a, double c) {
  size_t dim = 3;
  double center[] = {1.2, -3.4, 5.6};

  size_t num_directions = v2r_test_get_num_directions(dim);
  double *first = v2r_test_generate_directions(dim);
  double *last = first + dim * num_directions;

  double alpha_in = 0.9999;
  double alpha_out = 1.0001;
  char name[512];
  double ex[] = {1., 0., 0.}, p_in[dim], p_out[dim], d1[dim], d2[dim];
  for (double *d3 = first; d3 < last; d3 += dim) {
    V2R_Object *spheroid = v2r_spheroid_new(center, a, c, d3);
    v2r_cross(ex, d3, d1);
    v2r_normalize(d1);
    v2r_cross(d3, d1, d2);
    for (double *n = first; n < last; n += dim) {
      double n1 = v2r_dot(d1, n);
      double n2 = v2r_dot(d2, n);
      double n3 = v2r_dot(d3, n);
      for (size_t k = 0; k < dim; k++) {
        double r = a * (n1 * d1[k] + n2 * d2[k]) + c * n3 * d3[k];
        p_in[k] = spheroid->center[k] + alpha_in * r;
        p_out[k] = spheroid->center[k] + alpha_out * r;
      }
      sprintf(name, "/%s(a=%g,c=%g,d=(%g,%g,%g))/belongs/n=(%g,%g,%g)/true",
              spheroid->type->name, a, c, d3[0], d3[1], d3[2], n[0], n[1],
              n[2]);

      g_test_add_data_func_full(name,
                                v2r_test_belongs_data_new(spheroid, p_in, true),
                                v2r_test_belongs, v2r_test_belongs_data_free);

      sprintf(name, "/%s(a=%g,c=%g,d=(%g,%g,%g))/belongs/n=(%g,%g,%g)/false",
              spheroid->type->name, a, c, d3[0], d3[1], d3[2], n[0], n[1],
              n[2]);

      g_test_add_data_func_full(
          name, v2r_test_belongs_data_new(spheroid, p_out, false),
          v2r_test_belongs, v2r_test_belongs_data_free);
    }
    v2r_object_free(spheroid);
  }
  free(first);
}

void v2r_setup_test_spheroid() {
  g_test_add_func("/spheroid/new", v2r_test_spheroid_new);

  v2r_setup_test_spheroid_belongs(0.5, 0.02);
  v2r_setup_test_spheroid_belongs(0.02, 0.5);
}
