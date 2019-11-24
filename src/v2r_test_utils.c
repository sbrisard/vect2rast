#include <glib.h>

#include "v2r_test_utils.h"

void *v2r_test_belongs_data_new(V2R_Object *object, double const *point,
                                bool belongs) {
  V2R_TestBelongsData *data = malloc(sizeof(V2R_TestBelongsData));
  data->object = object;
  const size_t dim = data->object->type->dim;
  data->point = malloc(dim * sizeof(double));
  for (size_t i = 0; i < dim; i++) {
    data->point[i] = point[i];
  }
  data->belongs = belongs;
  return data;
}

void v2r_test_belongs_data_free(void *data_) {
  V2R_TestBelongsData *data = data_;
  v2r_object_free(data->object);
  free(data->point);
  free(data_);
}

void v2r_test_belongs(void const *data_) {
  V2R_TestBelongsData const *data = data_;
  g_assert_cmpint(data->object->type->belongs(data->object, data->point), ==, data->belongs);
}
