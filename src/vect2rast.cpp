#include <math.h>
#include <string.h>

#include "vect2rast/vect2rast.hpp"

static void init_bounds(double x_min, double x_max, double h_inv, int *i_min,
                        int *i_max) {
  *i_min = (int)ceil(h_inv * x_min - 0.5);
  *i_max = (int)floor(h_inv * x_max - 0.5);
}

///**
// * @brief Rasterize the specified object.
// *
// * For each cell of the grid whose center test_hypersphere_belongs to the object, the value is
// * set.
// *
// * @param object the object to rasterize
// * @param dim
// * @param size
// * @param grid
// * @param value
// */
//static int v2r_raster_3d(V2R_Object const *object, double const *length,
//                         size_t const *size, int *grid, int value) {
//  const size_t ndims = 3;
//  const double L0 = length[0], L1 = length[1], L2 = length[2];
//  // TODO This conversion is potentially unsagfe. However, we do need a signed
//  // value for periodic boundary conditions
//  const int n0 = size[0], n1 = size[1], n2 = size[2];
//  const double h0 = L0 / n0, h1 = L1 / n1, h2 = L2 / n2;
//
//  int i0_min, i0_max, i1_min, i1_max, i2_min, i2_max;
//  double bbmin[ndims], bbmax[ndims];
//  object->type->get_bounding_box(object, bbmin, bbmax);
//  init_bounds(bbmin[0], bbmax[0], n0 / L0, &i0_min, &i0_max);
//  init_bounds(bbmin[1], bbmax[1], n1 / L1, &i1_min, &i1_max);
//  init_bounds(bbmin[2], bbmax[2], n2 / L2, &i2_min, &i2_max);
//
//  double x[ndims];
//  int j0 = i0_min >= 0 ? i0_min : i0_min + n0;
//  x[0] = (i0_min + 0.5) * h0;
//  for (int i0 = i0_min; i0 <= i0_max; i0++) {
//    int j1 = i1_min >= 0 ? i1_min : i1_min + n1;
//    x[1] = (i1_min + 0.5) * h1;
//    for (int i1 = i1_min; i1 <= i1_max; i1++) {
//      int j2 = i2_min >= 0 ? i2_min : i2_min + n2;
//      x[2] = (i2_min + 0.5) * h2;
//      for (int i2 = i2_min; i2 <= i2_max; i2++) {
//        if (object->type->test_hypersphere_belongs(object, x)) {
//          grid[(j0 * n1 + j1) * n2 + j2] = value;
//        }
//        x[2] += h2;
//        ++j2;
//        if (j2 == n2) {
//          j2 = 0;
//        }
//      }
//      ++j1;
//      x[1] += h1;
//      if (j1 == n1) {
//        j1 = 0;
//      }
//    }
//    ++j0;
//    x[0] += h0;
//    if (j0 == n0) {
//      j0 = 0;
//    }
//  }
//  return 0;
//}
//
//static int v2r_raster_2d(V2R_Object const *object, double const *length,
//                         size_t const *size, int *grid, int value) {
//  const size_t ndims = 2;
//  const double L0 = length[0], L1 = length[1];
//  const int n0 = size[0], n1 = size[1];
//  const double h0 = L0 / n0, h1 = L1 / n1;
//
//  int i0_min, i0_max, i1_min, i1_max;
//  double bbmin[ndims], bbmax[ndims];
//  object->type->get_bounding_box(object, bbmin, bbmax);
//  init_bounds(bbmin[0], bbmax[0], n0 / L0, &i0_min, &i0_max);
//  init_bounds(bbmin[1], bbmax[1], n1 / L1, &i1_min, &i1_max);
//
//  double x[ndims];
//  int j0 = i0_min >= 0 ? i0_min : i0_min + n0;
//  x[0] = (i0_min + 0.5) * h0;
//  for (int i0 = i0_min; i0 <= i0_max; i0++) {
//    int j1 = i1_min >= 0 ? i1_min : i1_min + n1;
//    x[1] = (i1_min + 0.5) * h1;
//    for (int i1 = i1_min; i1 <= i1_max; i1++) {
//      if (object->type->test_hypersphere_belongs(object, x)) {
//        grid[j0 * n1 + j1] = value;
//      }
//      ++j1;
//      x[1] += h1;
//      if (j1 == n1) {
//        j1 = 0;
//      }
//    }
//    ++j0;
//    x[0] += h0;
//    if (j0 == n0) {
//      j0 = 0;
//    }
//  }
//  return 0;
//}
//
//int v2r_raster(V2R_Object const *object, double const *length,
//               size_t const *size, int *grid, int value) {
//  switch (object->type->dim) {
//    case 2:
//      return v2r_raster_2d(object, length, size, grid, value);
//    case 3:
//      return v2r_raster_3d(object, length, size, grid, value);
//    default:
//      return -1;
//  }
//}
