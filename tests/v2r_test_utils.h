#ifndef __V2R_TEST_UTILS_H_202001120826__
#define __V2R_TEST_UTILS_H_202001120826__
#include "vect2rast/vect2rast.h"

typedef struct V2R_TestRasterData {
  V2R_Object *object;
  double *length;
  size_t *size;
} V2R_TestRasterData;

DllExport V2R_TestRasterData *v2r_test_raster_data_new(V2R_Object *object,
                                                       double const *length,
                                                       size_t const *size);
DllExport void v2r_test_raster_data_free(void *data);
DllExport void v2r_test_raster(void const *data);

DllExport size_t v2r_test_get_num_directions(size_t dim);
DllExport double *v2r_test_generate_directions(size_t dim);

DllExport double v2r_dot(double const *v1, double const *v2);
DllExport void v2r_cross(double const *v1, double const *v2, double *v3);
DllExport void v2r_normalize(double *v);

DllExport void assert_equals_size_t(size_t expected, size_t actual);
DllExport void assert_equals_double(double expected, double actual, double rtol, double atol);

#endif
