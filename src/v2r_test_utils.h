#ifndef __V2R_TEST_UTILS_H__
#define __V2R_TEST_UTILS_H__
#include "vect2rast.h"

typedef struct V2R_TestBelongsData_ {
  V2R_Object *object;
  double *point;
  bool belongs;
} V2R_TestBelongsData;

DllExport void *v2r_test_belongs_data_new(V2R_Object const *object,
                                          double const *point, bool belongs);
DllExport void v2r_test_belongs_data_free(void *data);
DllExport void v2r_test_belongs(void const *data);

typedef struct V2R_TestRasterData {
  V2R_Object *object;
  double *length;
  size_t *size;
} V2R_TestRasterData;

DllExport V2R_TestRasterData *v2r_test_raster_data_new(V2R_Object *object,
                                                       double *const length,
                                                       size_t *const size);
DllExport void v2r_test_raster_data_free(void *data);
DllExport void v2r_test_raster(void const *data);

DllExport size_t v2r_test_get_num_directions(size_t dim);
DllExport double *v2r_test_generate_directions(size_t dim);

DllExport double v2r_dot(double const *v1, double const *v2);
DllExport void v2r_cross(double const *v1, double const *v2, double *v3);
DllExport void v2r_normalize(double *v);
#endif
