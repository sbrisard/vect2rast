#pragma once

#include "catch2/catch.hpp"

#include "test_utils.hpp"
#include "vect2rast/hypersphere.hpp"

template <size_t DIM>
void test_hypersphere_get_bounding_box(vect2rast::Hypersphere<DIM>& hypersphere) {
  std::array<double, DIM> bbmin, bbmax;
  hypersphere.get_bounding_box(bbmin, bbmax);

  for (size_t i = 0; i < DIM; i++) {
    REQUIRE(bbmin[i] == hypersphere.center[i] - hypersphere.radius);
    REQUIRE(bbmax[i] == hypersphere.center[i] + hypersphere.radius);
  }
}

template <size_t DIM>
void test_hypersphere_belongs(vect2rast::Hypersphere<DIM>& hypersphere) {
  auto directions = generate_directions<DIM>();
  double const r_in = 0.95 * hypersphere.radius;
  double const r_out = 1.05 * hypersphere.radius;

  std::array<double, DIM> p_in, p_out;
  for (const auto n : directions) {
    for (size_t j = 0; j < DIM; j++) {
      p_in[j] = hypersphere.center[j] + r_in * n[j];
      p_out[j] = hypersphere.center[j] + r_out * n[j];
    }

    REQUIRE(hypersphere.belongs(p_in));
    REQUIRE_FALSE(hypersphere.belongs(p_out));
  }
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

TEST_CASE("Hypersphere") {
  vect2rast::Hypersphere<2> disk{{1.2, -3.4}, 7.8};
  vect2rast::Hypersphere<3> sphere{{1.2, -3.4, 5.6}, 7.8};

  SECTION("Hypersphere.get_bounding_box") {
    test_hypersphere_get_bounding_box(disk);
    test_hypersphere_get_bounding_box(sphere);
  }

  SECTION("Hypersphere.belongs") {
    test_hypersphere_belongs(disk);
    test_hypersphere_belongs(sphere);
  }
  //  v2r_test_ndsphere_raster();
}
