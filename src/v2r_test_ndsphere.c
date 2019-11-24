#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "v2r_test_utils.h"
#include "vect2rast.h"

void test_disk_new() {
  size_t const dim = 2;
  double const c[] = {1.2, -3.4};
  double const r = 7.8;

  V2R_Object *disk = v2r_disk_new(c, r);
  g_assert_cmpint(disk->type->dim, ==, dim);
  g_assert_cmpfloat(v2r_ndsphere_radius(disk), ==, r);

  for (size_t i = 0; i < dim; i++) {
    g_assert_cmpfloat(disk->center[i], ==, c[i]);
    g_assert_cmpfloat(disk->bbmin[i], ==, c[i] - r);
    g_assert_cmpfloat(disk->bbmax[i], ==, c[i] + r);
  }

  v2r_object_free(disk);
}

void test_sphere_new() {
  size_t const dim = 3;
  double const c[] = {1.2, -3.4, 5.6};
  double const r = 7.8;

  V2R_Object *sphere = v2r_sphere_new(c, r);
  g_assert_cmpint(sphere->type->dim, ==, dim);
  g_assert_cmpfloat(v2r_ndsphere_radius(sphere), ==, r);

  for (size_t i = 0; i < dim; i++) {
    g_assert_cmpfloat(sphere->center[i], ==, c[i]);
    g_assert_cmpfloat(sphere->bbmin[i], ==, c[i] - r);
    g_assert_cmpfloat(sphere->bbmax[i], ==, c[i] + r);
  }

  v2r_object_free(sphere);
}

void v2r_setup_test_disk_belongs() {
  double const c[] = {1.2, -3.4};
  double const r = 7.8;
  V2R_Object *disk = v2r_disk_new(c, r);

  double *directions = v2r_test_generate_directions_2d();
  double *n = directions;
  double const r_in = 0.95 * r;
  double const r_out = 1.05 * r;

  char name[256];
  for (size_t i = 0; i < V2R_TEST_NUM_DIRECTIONS_2D; i++, n += 2) {
    double const p1[] = {c[0] + r_in * n[0], c[1] + r_in * n[1]};
    double const p2[] = {c[0] + r_out * n[0], c[1] + r_out * n[1]};

    sprintf(name, "/disk/belongs/in/%02d", (int)i);
    g_test_add_data_func_full(name, v2r_test_belongs_data_new(disk, p1, true),
                              v2r_test_belongs, v2r_test_belongs_data_free);

    sprintf(name, "/disk/belongs/out/%02d", (int)i);
    g_test_add_data_func_full(name, v2r_test_belongs_data_new(disk, p2, false),
                              v2r_test_belongs, v2r_test_belongs_data_free);
  }
  free(directions);
  v2r_object_free(disk);
}

void v2r_setup_test_sphere_belongs() {
  const size_t dim = 3;
  double const c[] = {1.2, -3.4, 5.6};
  double const r = 7.8;
  V2R_Object *sphere = v2r_sphere_new(c, r);
  double *directions = v2r_test_generate_directions_3d();
  double *n = directions;
  double const r_in = 0.95 * r;
  double const r_out = 1.05 * r;

  char name[256];
  double p_in[dim], p_out[dim];
  for (size_t i = 0; i < V2R_TEST_NUM_DIRECTIONS_3D; i++, n += dim) {
    for (size_t j = 0; j < dim; j++) {
      p_in[j] = c[j] + r_in * n[j];
      p_out[j] = c[j] + r_out * n[j];
    }

    sprintf(name, "/sphere/belongs/in/%02d", (int)i);
    g_test_add_data_func_full(name, v2r_test_belongs_data_new(sphere, p_in, true),
                              v2r_test_belongs, v2r_test_belongs_data_free);

    sprintf(name, "/sphere/belongs/out/%02d", (int)i);
    g_test_add_data_func_full(name,
                              v2r_test_belongs_data_new(sphere, p_out, false),
                              v2r_test_belongs, v2r_test_belongs_data_free);
  }
  v2r_object_free(sphere);
  free(directions);
}

void v2r_setup_test_ndsphere() {
  g_test_add_func("/disk/new", test_disk_new);
  g_test_add_func("/sphere/new", test_sphere_new);
  v2r_setup_test_disk_belongs();
  v2r_setup_test_sphere_belongs();
}
