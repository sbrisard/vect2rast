#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "v2r2.h"

void test_sphere_bbox(V2RObject const *sphere) {
  const double radius = V2R_OBJECT_DOUBLE_AT(sphere, V2R_SPHERE_RADIUS_INDEX);
  for (size_t i = 0; i < sphere->type->ndims; i++) {
    g_assert_cmpfloat(sphere->bbmin[i], ==, sphere->center[i] - radius);
    g_assert_cmpfloat(sphere->bbmax[i], ==, sphere->center[i] + radius);
  }
}

void setup_sphere_tests() {
  size_t ndims = 3;
  double c[] = {1.2, -3.4, 5.6};
  double r = 7.8;
  V2RObject const *sphere = v2r_sphere_new(ndims, c, r);

  V2RObject *sphere1 = v2r_object_copy(sphere);
  g_test_add_data_func_full("/sphere/bbox", sphere1, test_sphere_bbox,
                            v2r_object_free);
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  setup_sphere_tests();

  return g_test_run();
}
