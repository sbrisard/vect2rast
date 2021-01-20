#include <stdio.h>
#include <stdlib.h>

#include "v2r_test_utils.h"
#include "vect2rast/v2r_ndsphere.h"

void print_ndsphere(const V2R_Object *sphere) {
  printf("Sphere<%d>{r=%g, c=", sphere->type->dim,
         v2r_ndsphere_radius(sphere));
  print_array_double(sphere->type->dim, sphere->center);
  printf("}");
}

void test_disk_new() {
  printf("test_disk_new...");
  size_t const dim = 2;
  double const c[] = {1.2, -3.4};
  double const r = 7.8;

  V2R_Object *disk = v2r_disk_new(c, r);
  assert_equals_size_t(dim, disk->type->dim);
  assert_equals_double(r, v2r_ndsphere_radius(disk), 0.0, 0.0);

  for (size_t i = 0; i < dim; i++) {
    assert_equals_double(c[i], disk->center[i], 0.0, 0.0);
  }

  v2r_object_free(disk);
  printf(" OK\n");
}

void test_disk_get_bounding_box() {
  printf("test_disk_bounding_box...");
  size_t const dim = 2;
  double const c[] = {1.2, -3.4};
  double const r = 7.8;

  V2R_Object *disk = v2r_disk_new(c, r);

  double bbmin[dim], bbmax[dim];
  disk->type->get_bounding_box(disk, bbmin, bbmax);

  for (size_t i = 0; i < dim; i++) {
    assert_equals_double(c[i] - r, bbmin[i], 0.0, 0.0);
    assert_equals_double(c[i] + r, bbmax[i], 0.0, 0.0);
  }

  v2r_object_free(disk);
  printf(" OK\n");
}

void test_sphere_new() {
  printf("test_sphere_new...");
  size_t const dim = 3;
  double const c[] = {1.2, -3.4, 5.6};
  double const r = 7.8;

  V2R_Object *sphere = v2r_sphere_new(c, r);
  assert_equals_size_t(dim, sphere->type->dim);
  assert_equals_double(r, v2r_ndsphere_radius(sphere), 0.0, 0.0);

  for (size_t i = 0; i < dim; i++) {
    assert_equals_double(c[i], sphere->center[i], 0.0, 0.0);
  }

  v2r_object_free(sphere);
  printf(" OK\n");
}

void test_sphere_get_bounding_box() {
  printf("test_sphere_get_bounding_box...");
  size_t const dim = 3;
  double const c[] = {1.2, -3.4, 5.6};
  double const r = 7.8;

  V2R_Object *sphere = v2r_sphere_new(c, r);

  double bbmin[dim], bbmax[dim];
  sphere->type->get_bounding_box(sphere, bbmin, bbmax);

  for (size_t i = 0; i < dim; i++) {
    assert_equals_double(c[i] - r, bbmin[i], 0.0, 0.0);
    assert_equals_double(c[i] + r, bbmax[i], 0.0, 0.0);
  }

  v2r_object_free(sphere);
  printf(" OK\n");
}

void v2r_test_ndsphere_belongs(V2R_Object const *ndsphere) {
  printf("test_ndsphere_belongs ");
  print_ndsphere(ndsphere);
  printf("...");
  const size_t dim = ndsphere->type->dim;
  double *directions = v2r_test_generate_directions(dim);
  double const *n = directions;
  double const r = v2r_ndsphere_radius(ndsphere);
  double const r_in = 0.95 * r;
  double const r_out = 1.05 * r;

  double p_in[dim], p_out[dim];
  for (size_t i = 0; i < v2r_test_get_num_directions(dim); i++, n += dim) {
    for (size_t j = 0; j < dim; j++) {
      p_in[j] = ndsphere->center[j] + r_in * n[j];
      p_out[j] = ndsphere->center[j] + r_out * n[j];
    }

    assert_true(ndsphere->type->belongs(ndsphere, p_in));
    assert_false(ndsphere->type->belongs(ndsphere, p_out));
  }
  free(directions);
  printf(" OK\n");
}

void v2r_test_ndsphere_raster() {
  double length[] = {1.5, 2.6, 3.7};
  size_t size[] = {50, 60, 70};
  double xi = 0.05;
  double eta = 0.06;
  double zeta = 0.07;
  double c[3];
  double r = 0.5;

  for (size_t i0 = 0; i0 <= 1; i0++) {
    c[0] = i0 == 0 ? xi * length[0] : (1. - xi) * length[0];
    for (size_t i1 = 0; i1 <= 1; i1++) {
      c[1] = i1 == 0 ? eta * length[1] : (1. - eta) * length[1];
      for (size_t i2 = 0; i2 <= 1; i2++) {
        c[2] = i2 == 0 ? zeta * length[2] : (1. - zeta) * length[2];
        V2R_Object *sphere = v2r_sphere_new(c, r);
        v2r_test_raster(sphere, length, size);
        v2r_object_free(sphere);
      }
    }
  }
}

void v2r_test_ndsphere_all() {
  test_disk_new();
  test_disk_get_bounding_box();
  test_sphere_new();
  test_sphere_get_bounding_box();

  double const c[] = {1.2, -3.4, 5.6};
  double const r = 7.8;
  V2R_Object *disk = v2r_disk_new(c, r);
  v2r_test_ndsphere_belongs(disk);
  v2r_object_free(disk);
  V2R_Object *sphere = v2r_sphere_new(c, r);
  v2r_test_ndsphere_belongs(sphere);
  v2r_object_free(sphere);

  v2r_test_ndsphere_raster();
}
