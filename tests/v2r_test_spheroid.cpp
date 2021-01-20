#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "v2r_test_utils.hpp"
#include "vect2rast/v2r_spheroid.hpp"

void v2r_test_spheroid_new() {
  printf("test_spheroid_new...");
  size_t const dim = 3;
  double const x[] = {1.2, -3.4, 5.6};
  double const theta = .3 * M_PI;
  double const phi = .5 * M_PI;
  double const n[] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
  double const a = 7.8;
  double const c = 9.1;

  V2R_Object *spheroid = v2r_spheroid_new(x, a, c, n);
  assert_equals_size_t(dim, spheroid->type->dim);
  assert_equals_double(a, v2r_spheroid_equatorial_radius(spheroid), 0., 0.);
  assert_equals_double(c, v2r_spheroid_polar_radius(spheroid), 0., 0.);

  double n_act[3];
  v2r_spheroid_axis(spheroid, n_act);
  for (size_t i = 0; i < dim; i++) {
    assert_equals_double(x[i], spheroid->center[i], 0., 0.);
    assert_equals_double(n[i], n_act[i], 0., 0.);
  }

  v2r_object_free(spheroid);
  printf(" OK\n");
}

void v2r_test_spheroid_belongs(double *const center, double a, double c) {
  size_t const dim = 3;
  printf("test_spheroid_belongs{center=");
  print_array_double(dim, center);
  printf(", a=%g, c=%g}... ", a, c);

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
      assert_true(spheroid->type->belongs(spheroid, p_in));
      assert_false(spheroid->type->belongs(spheroid, p_out));
    }
    v2r_object_free(spheroid);
  }
  free(first);
  printf("OK\n");
}

void v2r_test_spheroid_all() {
  v2r_test_spheroid_new();

  double center[] = {1.2, -3.4, 5.6};
  v2r_test_spheroid_belongs(center, 0.5, 0.02);
  v2r_test_spheroid_belongs(center, 0.02, 0.5);
}
