#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "v2r_test_utils.hpp"
#include "vect2rast/v2r_ndsphere.hpp"

template <size_t DIM>
void test_hypersphere_get_bounding_box(Hypersphere<DIM> hypersphere) {
  std::cout << "test_hypersphere_get_bounding_box(" << hypersphere << ")...";

  double bbmin[DIM], bbmax[DIM];
  hypersphere.get_bounding_box(bbmin, bbmax);

  for (size_t i = 0; i < DIM; i++) {
    assert_equals_double(hypersphere.center[i] - hypersphere.radius, bbmin[i],
                         0.0, 0.0);
    assert_equals_double(hypersphere.center[i] + hypersphere.radius, bbmax[i],
                         0.0, 0.0);
  }

  std::cout << "OK" << std::endl;
}

template <size_t DIM>
void test_hypersphere_belongs(Hypersphere<DIM> hypersphere) {
  std::cout << "test_hypersphere_belongs(" << hypersphere << ")...";
  double *directions = v2r_test_generate_directions(DIM);
  double const *n = directions;
  double const r_in = 0.95 * hypersphere.radius;
  double const r_out = 1.05 * hypersphere.radius;

  double p_in[DIM], p_out[DIM];
  for (size_t i = 0; i < v2r_test_get_num_directions(DIM); i++, n += DIM) {
    for (size_t j = 0; j < DIM; j++) {
      p_in[j] = hypersphere.center[j] + r_in * n[j];
      p_out[j] = hypersphere.center[j] + r_out * n[j];
    }

    assert_true(hypersphere.belongs(p_in));
    assert_false(hypersphere.belongs(p_out));
  }
  free(directions);
  std::cout << "OK" << std::endl;
}

// void v2r_test_ndsphere_raster() {
//  double length[] = {1.5, 2.6, 3.7};
//  size_t size[] = {50, 60, 70};
//  double xi = 0.05;
//  double eta = 0.06;
//  double zeta = 0.07;
//  double c[3];
//  double r = 0.5;
//
//  for (size_t i0 = 0; i0 <= 1; i0++) {
//    c[0] = i0 == 0 ? xi * length[0] : (1. - xi) * length[0];
//    for (size_t i1 = 0; i1 <= 1; i1++) {
//      c[1] = i1 == 0 ? eta * length[1] : (1. - eta) * length[1];
//      for (size_t i2 = 0; i2 <= 1; i2++) {
//        c[2] = i2 == 0 ? zeta * length[2] : (1. - zeta) * length[2];
//        V2R_Object *sphere = v2r_sphere_new(c, r);
//        v2r_test_raster(sphere, length, size);
//        v2r_object_free(sphere);
//      }
//    }
//  }
//}

void test_hypersphere_all() {
  Hypersphere<2> disk{{1.2, -3.4}, 7.8};
  Hypersphere<3> sphere{{1.2, -3.4, 5.6}, 7.8};

  test_hypersphere_get_bounding_box(disk);
  test_hypersphere_get_bounding_box(sphere);

  test_hypersphere_belongs(disk);
  test_hypersphere_belongs(sphere);

  //  v2r_test_ndsphere_raster();
}
