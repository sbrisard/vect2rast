#include <glib.h>
#include <math.h>

#include "v2r_test_utils.h"

void *v2r_test_belongs_data_new(V2R_Object const *object, double const *point,
                                bool belongs) {
  V2R_TestBelongsData *data = malloc(sizeof(V2R_TestBelongsData));
  data->object = v2r_object_copy(object);
  const size_t size = object->type->dim * sizeof(double);
  data->point = malloc(size);
  memcpy(data->point, point, size);
  data->belongs = belongs;
  return data;
}

void v2r_test_belongs_data_free(void *data) {
  V2R_TestBelongsData *data_ = data;
  v2r_object_free(data_->object);
  free(data_->point);
  free(data);
}

void v2r_test_belongs(void const *data) {
  V2R_TestBelongsData const *data_ = data;
  g_assert_cmpuint(data_->object->type->belongs(data_->object, data_->point),
                   ==, data_->belongs);
}

V2R_TestRasterData *v2r_test_raster_data_new(V2R_Object *object,
                                             double *const length,
                                             size_t *const size) {
  size_t const dim = object->type->dim;
  V2R_TestRasterData *data = malloc(sizeof(V2R_TestRasterData));
  data->object = v2r_object_copy(object);
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
