#include <glib.h>
#include <math.h>

#include "v2r_test_utils.h"

void *v2r_test_belongs_data_new(V2R_Object *object, double const *point,
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

double *v2r_test_generate_directions_2d() {
  double *directions = malloc(2 * V2R_TEST_NUM_DIRECTIONS_2D * sizeof(double));
  double *dir = directions;
  for (size_t i = 0; i < V2R_TEST_NUM_DIRECTIONS_2D; i++) {
    double const theta = 2 * M_PI * i / (double)V2R_TEST_NUM_DIRECTIONS_2D;
    *dir = cos(theta);
    dir += 1;
    *dir = sin(theta);
    dir += 1;
  }
  return directions;
}
