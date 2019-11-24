#include <glib.h>

#include "v2r_test_utils.h"

void *v2r_test_belongs_data_new(V2R_Object *object, double const *point,
                                bool belongs) {
  V2R_TestBelongsData *data = malloc(sizeof(V2R_TestBelongsData));
  data->object = object;
  data->point = point;
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
  g_assert_cmpuint(data_->object->type->belongs(data_->object, data_->point), ==,
                    data_->belongs);
}
