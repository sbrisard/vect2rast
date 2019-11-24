#ifndef __V2R_TEST_UTILS_H__
#define __V2R_TEST_UTILS_H__
#include "vect2rast.h"

#define V2R_TEST_NUM_DIRECTIONS_2D 10

typedef struct V2R_TestBelongsData_ {
  V2R_Object *object;
  double *point;
  bool belongs;
} V2R_TestBelongsData;

DllExport void *v2r_test_belongs_data_new(V2R_Object *object,
                                          double const *point, bool belongs);
DllExport void v2r_test_belongs_data_free(void *data);
DllExport void v2r_test_belongs(void const *data);

DllExport double *v2r_test_generate_directions_2d();

#endif
