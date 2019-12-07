#include <math.h>
#include <string.h>

#include "vect2rast.h"

V2R_Object *v2r_object_new(V2R_ObjectType const *type) {
  size_t const dim = type->dim;
  V2R_Object *object = malloc(sizeof(V2R_Object));
  object->type = type;

  size_t const size = dim * sizeof(double);
  object->center = malloc(size);
  object->bbmin = malloc(size);
  object->bbmax = malloc(size);

  object->data = NULL;

  return object;
}

void v2r_object_free(V2R_Object *object) {
  object->type->data_free(object->data);
  free(object->center);
  free(object->bbmin);
  free(object->bbmax);
}

V2R_Object *v2r_object_copy(V2R_Object const *object, double const *center) {
  size_t const size = object->type->dim * sizeof(double);
  V2R_Object *copy = v2r_object_new(object->type);
  copy->data = object->type->data_copy(object->data);
  if (center == NULL) {
    memcpy(copy->center, object->center, size);
    memcpy(copy->bbmin, object->bbmin, size);
    memcpy(copy->bbmax, object->bbmax, size);
  } else {
    memcpy(copy->center, center, size);
    copy->type->init_bounding_box(copy);
  }
  return copy;
}

static void init_bounds(double x_min, double x_max, double h_inv, int *i_min,
                        int *i_max) {
  *i_min = (int)ceil(h_inv * x_min - 0.5);
  *i_max = (int)floor(h_inv * x_max - 0.5);
}

/**
 * @brief Rasterize the specified object.
 *
 * For each cell of the grid whose center belongs to the object, the value is
 * set.
 *
 * @param object the object to rasterize
 * @param dim
 * @param size
 * @param grid
 * @param value
 */
static int v2r_raster_3d(V2R_Object *object, double const *length,
                         size_t const *size, int *grid, int value) {
  const size_t ndims = 3;
  const double L0 = length[0], L1 = length[1], L2 = length[2];
  const int n0 = size[0], n1 = size[1], n2 = size[2];
  const double h0 = L0 / n0, h1 = L1 / n1, h2 = L2 / n2;

  int i0_min, i0_max, i1_min, i1_max, i2_min, i2_max;
  init_bounds(object->bbmin[0], object->bbmax[0], n0 / L0, &i0_min, &i0_max);
  init_bounds(object->bbmin[1], object->bbmax[1], n1 / L1, &i1_min, &i1_max);
  init_bounds(object->bbmin[2], object->bbmax[2], n2 / L2, &i2_min, &i2_max);

  double x[ndims];
  int j0 = i0_min >= 0 ? i0_min : i0_min + n0;
  x[0] = (i0_min + 0.5) * h0;
  for (int i0 = i0_min; i0 <= i0_max; i0++) {
    int j1 = i1_min >= 0 ? i1_min : i1_min + n1;
    x[1] = (i1_min + 0.5) * h1;
    for (int i1 = i1_min; i1 <= i1_max; i1++) {
      int j2 = i2_min >= 0 ? i2_min : i2_min + n2;
      x[2] = (i2_min + 0.5) * h2;
      for (int i2 = i2_min; i2 <= i2_max; i2++) {
        if (object->type->belongs(object, x)) {
          grid[(j0 * n1 + j1) * n2 + j2] = value;
        }
        x[2] += h2;
        ++j2;
        if (j2 == n2) {
          j2 = 0;
        }
      }
      ++j1;
      x[1] += h1;
      if (j1 == n1) {
        j1 = 0;
      }
    }
    ++j0;
    x[0] += h0;
    if (j0 == n0) {
      j0 = 0;
    }
  }
  return 0;
}

static int v2r_raster_2d(V2R_Object *object, double const *length,
                         size_t const *size, int *grid, int value) {
  const size_t ndims = 2;
  const double L0 = length[0], L1 = length[1];
  const int n0 = size[0], n1 = size[1];
  const double h0 = L0 / n0, h1 = L1 / n1;

  int i0_min, i0_max, i1_min, i1_max;
  init_bounds(object->bbmin[0], object->bbmax[0], n0 / L0, &i0_min, &i0_max);
  init_bounds(object->bbmin[1], object->bbmax[1], n1 / L1, &i1_min, &i1_max);

  double x[ndims];
  int j0 = i0_min >= 0 ? i0_min : i0_min + n0;
  x[0] = (i0_min + 0.5) * h0;
  for (int i0 = i0_min; i0 <= i0_max; i0++) {
    int j1 = i1_min >= 0 ? i1_min : i1_min + n1;
    x[1] = (i1_min + 0.5) * h1;
    for (int i1 = i1_min; i1 <= i1_max; i1++) {
      if (object->type->belongs(object, x)) {
        grid[j0 * n1 + j1] = value;
      }
      ++j1;
      x[1] += h1;
      if (j1 == n1) {
        j1 = 0;
      }
    }
    ++j0;
    x[0] += h0;
    if (j0 == n0) {
      j0 = 0;
    }
  }
  return 0;
}

int v2r_raster(V2R_Object *object, double *length, size_t *size, int *grid,
               int value) {
  switch (object->type->dim) {
    case 2:
      return v2r_raster_2d(object, length, size, grid, value);
    case 3:
      return v2r_raster_3d(object, length, size, grid, value);
    default:
      return -1;
  }
}
