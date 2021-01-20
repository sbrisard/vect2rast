#ifndef __V2R_TEST_UTILS_H_202001120826__
#define __V2R_TEST_UTILS_H_202001120826__
#include <stdio.h>
#include <string.h>
#include "vect2rast/vect2rast.h"

DllExport void v2r_test_raster(V2R_Object const *object, double const *length,
                               size_t const *size);

DllExport size_t v2r_test_get_num_directions(size_t dim);
DllExport double *v2r_test_generate_directions(size_t dim);

DllExport double v2r_dot(double const *v1, double const *v2);
DllExport void v2r_cross(double const *v1, double const *v2, double *v3);
DllExport void v2r_normalize(double *v);

DllExport void print_array_size_t(size_t n, const size_t *a);
DllExport void print_array_double(size_t n, const double *a);

DllExport void assert_true(bool predicate);
DllExport void assert_false(bool predicate);
DllExport void assert_equals_size_t(size_t expected, size_t actual);
DllExport void assert_equals_double(double expected, double actual, double rtol,
                                    double atol);

#endif
