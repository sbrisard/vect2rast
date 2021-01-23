#include "test_utils.hpp"
#include <iostream>

void assert_true(bool predicate) {
  if (!predicate) exit(-1);
}

void assert_false(bool predicate) {
  if (predicate) exit(-1);
}

void assert_equals_double(double expected, double actual, double rtol,
                          double atol) {
  if (fabs(actual - expected) > rtol * fabs(expected) + atol) {
    exit(-1);
  }
}

void minimum_image(double L, double L_half, double* x) {
  if (*x < -L_half) *x += L;
  if (*x > L_half) *x -= L;
}

// void v2r_test_raster(V2R_Object const *object, double const *length,
//                     size_t const *size) {
//  std::cout << "v2r_test_raster{" << object->type->name << "}<"
//            << object->type->dim << ">";
//  print_array_double(object->type->dim, object->center);
//  printf("}, size=");
//  print_array_size_t(object->type->dim, size);
//  printf(", length=");
//  print_array_double(object->type->dim, length);
//  printf("...");
//  size_t const max_dim = 3;
//  size_t const dim = object->type->dim;
//
//  size_t n[] = {size[0], size[1], 1};
//  double L[] = {length[0], length[1], 0.};
//  double c[] = {object->center[0], object->center[1], 0.};
//  if (dim == max_dim) {
//    size_t const i = dim - 1;
//    n[i] = size[i];
//    L[i] = length[i];
//    c[i] = object->center[i];
//  }
//  double h[max_dim], L_half[max_dim];
//  for (size_t i = 0; i < max_dim; i++) {
//    h[i] = L[i] / n[i];
//    L_half[i] = 0.5 * L[i];
//  }
//
//  auto actual = static_cast<int *>(calloc(n[0] * n[1] * n[2], sizeof(int)));
//  v2r_raster(object, length, size, actual, 1);
//
//  /* Create copy of particle, centered at the origin, since we will compute
//  the
//   * minimum image of the CP vector (C: center; P: current point). */
//  double O[] = {0., 0., 0.};
//  V2R_Object *object2 = v2r_object_copy(object, O);
//  double point[dim];
//  for (size_t i0 = 0; i0 < n[0]; i0++) {
//    point[0] = (i0 + 0.5) * h[0] - c[0];
//    minimum_image(L[0], L_half[0], point);
//    for (size_t i1 = 0; i1 < n[1]; i1++) {
//      point[1] = (i1 + 0.5) * h[1] - c[1];
//      minimum_image(L[1], L_half[1], point + 1);
//      for (size_t i2 = 0; i2 < n[2]; i2++) {
//        point[2] = (i2 + 0.5) * h[2] - c[2];
//        minimum_image(L[2], L_half[2], point + 2);
//        size_t j = (i0 * n[1] + i1) * n[2] + i2;
//        int expected = object2->type->belongs(object2, point);
//        assert_equals_int(expected, actual[j]);
//      }
//    }
//  }
//  free(actual);
//  printf(" OK\n");
//}

size_t v2r_test_get_num_directions(size_t dim) {
  switch (dim) {
    case 2:
      return 10;
    case 3:
      return 12;
    default:
      return 0;
  }
}

std::array<double, 3> v2r_cross(const std::array<double, 3>& v1,
                                const std::array<double, 3>& v2) {
  return {v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2],
          v1[0] * v2[1] - v1[1] * v2[0]};
}

void v2r_normalize(std::array<double, 3>& v) {
  double s = 1. / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}
