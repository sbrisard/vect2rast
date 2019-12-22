#include <glib.h>
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
  }

  v2r_object_free(disk);
}

void test_disk_get_bounding_box() {
  size_t const dim = 2;
  double const c[] = {1.2, -3.4};
  double const r = 7.8;

  V2R_Object *disk = v2r_disk_new(c, r);

  double bbmin[dim], bbmax[dim];
  disk->type->get_bounding_box(disk, bbmin, bbmax);

  for (size_t i = 0; i < dim; i++) {
    g_assert_cmpfloat(bbmin[i], ==, c[i] - r);
    g_assert_cmpfloat(bbmax[i], ==, c[i] + r);
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
  }

  v2r_object_free(sphere);
}

void test_sphere_get_bounding_box() {
  size_t const dim = 3;
  double const c[] = {1.2, -3.4, 5.6};
  double const r = 7.8;

  V2R_Object *sphere = v2r_sphere_new(c, r);

  double bbmin[dim], bbmax[dim];
  sphere->type->get_bounding_box(sphere, bbmin, bbmax);

  for (size_t i = 0; i < dim; i++) {
    g_assert_cmpfloat(bbmin[i], ==, c[i] - r);
    g_assert_cmpfloat(bbmax[i], ==, c[i] + r);
  }

  v2r_object_free(sphere);
}

void v2r_setup_test_ndsphere_belongs(V2R_Object const *ndsphere) {
  const size_t dim = ndsphere->type->dim;
  double *directions = v2r_test_generate_directions(dim);
  double const *n = directions;
  double const r = v2r_ndsphere_radius(ndsphere);
  double const r_in = 0.95 * r;
  double const r_out = 1.05 * r;

  char name[256];
  double p_in[dim], p_out[dim];
  for (size_t i = 0; i < v2r_test_get_num_directions(dim); i++, n += dim) {
    for (size_t j = 0; j < dim; j++) {
      p_in[j] = ndsphere->center[j] + r_in * n[j];
      p_out[j] = ndsphere->center[j] + r_out * n[j];
    }

    sprintf(name, "/%s/belongs/in/%02d", ndsphere->type->name, (int) i);
    g_test_add_data_func_full(name,
                              v2r_test_belongs_data_new(ndsphere, p_in, true),
                              v2r_test_belongs, v2r_test_belongs_data_free);

    sprintf(name, "/%s/belongs/out/%02d", ndsphere->type->name, (int) i);
    g_test_add_data_func_full(name,
                              v2r_test_belongs_data_new(ndsphere, p_out, false),
                              v2r_test_belongs, v2r_test_belongs_data_free);
  }
  free(directions);
}

void v2r_setup_test_ndsphere_raster() {
  double length[] = {1.5, 2.6, 3.7};
  size_t size[] = {50, 60, 70};
  double xi = 0.05;
  double eta = 0.06;
  double zeta = 0.07;
  double c[3];
  double r = 0.5;
  char name[255];

  for (size_t i0 = 0; i0 <= 1; i0++) {
    c[0] = i0 == 0 ? xi * length[0] : (1. - xi) * length[0];
    for (size_t i1 = 0; i1 <= 1; i1++) {
      c[1] = i1 == 0 ? eta * length[1] : (1. - eta) * length[1];
      for (size_t i2 = 0; i2 <= 1; i2++) {
        c[2] = i2 == 0 ? zeta * length[2] : (1. - zeta) * length[2];
        sprintf(name, "/Sphere/raster/%zu", (i0 * 2 + i1) * 2 + i2);
        V2R_Object *sphere = v2r_sphere_new(c, r);
        g_test_add_data_func_full(
            name, v2r_test_raster_data_new(sphere, length, size),
            v2r_test_raster, v2r_test_raster_data_free);
        v2r_object_free(sphere);
      }
    }
  }
}

void v2r_setup_test_ndsphere() {
  g_test_add_func("/Disk/new", test_disk_new);
  g_test_add_func("/Disk/get_bounding_box", test_disk_get_bounding_box);
  g_test_add_func("/Sphere/new", test_sphere_new);
  g_test_add_func("/Sphere/get_bounding_box", test_sphere_get_bounding_box);

  double const c[] = {1.2, -3.4, 5.6};
  double const r = 7.8;
  V2R_Object *disk = v2r_disk_new(c, r);
  V2R_Object *sphere = v2r_sphere_new(c, r);
  v2r_setup_test_ndsphere_belongs(disk);
  v2r_setup_test_ndsphere_belongs(sphere);
  v2r_setup_test_ndsphere_raster();
  v2r_object_free(disk);
  v2r_object_free(sphere);
}
