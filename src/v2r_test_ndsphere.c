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

void test_disk_belongs() {
  double const c[] = {1.2, -3.4};
  double const r = 7.8;
  V2R_Object *disk = v2r_disk_new(c, r);

  size_t const num_points = 10;

  double const r1 = 0.95 * r;
  double const r2 = 1.05 * r;

  for (size_t i = 0; i < num_points; i++) {
    double const theta = 2 * M_PI * i / (double)num_points;
    double const sin_theta = sin(theta);
    double const cos_theta = cos(theta);
    double const p1[] = {c[0] + r1 * cos_theta, c[1] + r1 * sin_theta};
    g_assert_cmpuint(disk->type->belongs(disk, p1), ==, true);

    double const p2[] = {c[0] + r2 * cos_theta, c[1] + r2 * sin_theta};
    g_assert_cmpuint(disk->type->belongs(disk, p2), ==, false);
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

void test_sphere_belongs() {
  double const c[] = {1.2, -3.4, 5.6};
  double const r = 7.8;
  V2R_Object *sphere = v2r_sphere_new(c, r);

  size_t const num_points = 10;

  double const r1 = 0.95 * r;
  double const r2 = 1.05 * r;

  for (size_t i = 0; i < num_points; i++) {
    double const theta = M_PI * i / (num_points - 1.);
    double const sin_theta = sin(theta);
    double const cos_theta = cos(theta);
    for (size_t j = 0; j < num_points; j++) {
      double const phi = 2 * M_PI * i / (double)num_points;
      double const sin_phi = sin(phi);
      double const cos_phi = cos(phi);

      double const nx = sin_theta * cos_phi;
      double const ny = sin_theta * sin_phi;
      double const nz = cos_theta;

      double const p1[] = {c[0] + r1 * nx, c[1] + r1 * ny, c[2] + r1 * nz};
      g_assert_cmpuint(sphere->type->belongs(sphere, p1), ==, true);

      double const p2[] = {c[0] + r2 * nx, c[1] + r2 * ny, c[2] + r2 * nz};
      g_assert_cmpuint(sphere->type->belongs(sphere, p2), ==, false);
    }
  }

  v2r_object_free(sphere);
}

void v2r_setup_test_disk_belongs() {
  double const c[] = {1.2, -3.4};
  double const r = 7.8;
  V2R_Object *disk = v2r_disk_new(c, r);

  size_t const num_points = 10;

  double const r_in = 0.95 * r;
  double const r_out = 1.05 * r;

  char name[256];
  for (size_t i = 0; i < num_points; i++) {
    double const theta = 2 * M_PI * i / (double)num_points;
    double const sin_theta = sin(theta);
    double const cos_theta = cos(theta);
    double const p1[] = {c[0] + r_in * cos_theta, c[1] + r_in * sin_theta};
    double const p2[] = {c[0] + r_out * cos_theta, c[1] + r_out * sin_theta};

    sprintf(name, "/disk/belongs/in/%02d", (int)i);
    g_test_add_data_func_full(name, v2r_test_belongs_data_new(disk, p1, true),
                              v2r_test_belongs, v2r_test_belongs_data_free);

    sprintf(name, "/disk/belongs/out/%02d", (int)i);
    g_test_add_data_func_full(name, v2r_test_belongs_data_new(disk, p2, false),
                              v2r_test_belongs, v2r_test_belongs_data_free);
  }
  v2r_object_free(disk);
}

void v2r_setup_test_ndsphere() {
  g_test_add_func("/disk/new", test_disk_new);
  g_test_add_func("/sphere/new", test_sphere_new);
  g_test_add_func("/sphere/belongs", test_sphere_belongs);
  v2r_setup_test_disk_belongs();
}
