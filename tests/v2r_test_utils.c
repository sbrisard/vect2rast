#include <glib.h>
#include <math.h>

#include "v2r_test_utils.h"

V2R_TestRasterData *v2r_test_raster_data_new(V2R_Object *object,
                                             double const *length,
                                             size_t const *size) {
  size_t const dim = object->type->dim;
  V2R_TestRasterData *data = malloc(sizeof(V2R_TestRasterData));
  data->object = v2r_object_copy(object, NULL);
  data->length = malloc(dim * sizeof(double));
  data->size = malloc(dim * sizeof(size_t));
  for (size_t i = 0; i < dim; i++) {
    data->length[i] = length[i];
    data->size[i] = size[i];
  }
  return data;
}

void v2r_test_raster_data_free(void *data) {
  V2R_TestRasterData *data_ = data;
  v2r_object_free(data_->object);
  free(data_->size);
  free(data_->length);
}

void minimum_image(double L, double L_half, double *x) {
  if (*x < -L_half)
    *x += L;
  if (*x > L_half)
    *x -= L;
}

void v2r_test_raster(void const * data) {
  V2R_TestRasterData const *data_ = data;
  size_t const max_dim = 3;
  size_t const dim = data_->object->type->dim;

  size_t n[] = {data_->size[0], data_->size[1], 1};
  double L[] = {data_->length[0], data_->length[1], 0.};
  double c[] = {data_->object->center[0], data_->object->center[1], 0.};
  if (dim == max_dim) {
    size_t const i = dim - 1;
    n[i] = data_->size[i];
    L[i] = data_->length[i];
    c[i] = data_->object->center[i];
  }
  double h[max_dim], L_half[max_dim];
  for (size_t i = 0; i < max_dim; i++) {
    h[i] = L[i] / n[i];
    L_half[i] = 0.5 * L[i];
  }

  int *actual = calloc(n[0] * n[1] * n[2], sizeof(int));
  v2r_raster(data_->object, data_->length, data_->size, actual, 1);

  /* Create copy of particle, centered at the origin, since we will compute the
   * minimum image of the CP vector (C: center; P: current point). */
  double O[] = {0., 0., 0.};
  V2R_Object *object = v2r_object_copy(data_->object, O);
  double point[dim];
  for (size_t i0 = 0; i0 < n[0]; i0++) {
    point[0] = (i0 + 0.5) * h[0] - c[0];
    minimum_image(L[0], L_half[0], point);
    for (size_t i1 = 0; i1 < n[1]; i1++) {
      point[1] = (i1 + 0.5) * h[1] - c[1];
      minimum_image(L[1], L_half[1], point + 1);
      for (size_t i2 = 0; i2 < n[2]; i2++) {
        point[2] = (i2 + 0.5) * h[2] - c[2];
        minimum_image(L[2], L_half[2], point + 2);
        size_t j = (i0 * n[1] + i1) * n[2] + i2;
        int expected = object->type->belongs(object, point);
        g_assert_cmpint(expected, ==, actual[j]);
      }
    }
  }
  free(actual);
}

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

double *v2r_test_generate_directions(size_t dim) {
  size_t const num_directions = v2r_test_get_num_directions(dim);
  size_t const size = dim * num_directions * sizeof(double);
  double *directions = malloc(size);
  if (dim == 2) {
    double *dir = directions;
    for (size_t i = 0; i < num_directions; i++) {
      double const theta = 2 * M_PI * i / (double)num_directions;
      *dir = cos(theta);
      dir += 1;
      *dir = sin(theta);
      dir += 1;
    }
  } else if (dim == 3) {
    double phi = .5 * (1. + sqrt(5.));
    double u = 1. / sqrt(1 + phi * phi);
    double v = phi * u;
    double dir[] = {0., -u, -v, 0., -u, +v, 0., +u, -v, 0., +u, +v,
                    -u, -v, 0., -u, +v, 0., +u, -v, 0., +u, +v, 0.,
                    -v, 0., -u, -v, 0., +u, +v, 0., -u, +v, 0., +u};
    memcpy(directions, dir, size);
  } else {
    return NULL;
  }
  return directions;
}

double v2r_dot(double const *v1, double const *v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void v2r_cross(double const *v1, double const *v2, double *v3) {
  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void v2r_normalize(double *v) {
  double s = 1. / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}
