#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <vect2rast.h>

void test_ndsphere_bbox(V2R_Object const *sphere) {
  const double radius = v2r_ndsphere_radius(sphere);
  for (size_t i = 0; i < sphere->type->ndims; i++) {
    g_assert_cmpfloat(sphere->bbmin[i], ==, sphere->center[i] - radius);
    g_assert_cmpfloat(sphere->bbmax[i], ==, sphere->center[i] + radius);
  }
}

void test_disk_setup_tests() {
  double c[] = {1.2, -3.4};
  double r = 7.8;

  g_test_add_data_func_full("/disk/bbox", v2r_disk_new(c, r),
                            test_ndsphere_bbox, v2r_object_free);
}

void test_sphere_setup_tests() {
  double c[] = {1.2, -3.4, 5.6};
  double r = 7.8;

  g_test_add_data_func_full("/sphere/bbox", v2r_sphere_new(c, r),
                            test_ndsphere_bbox, v2r_object_free);
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  test_disk_setup_tests();
  test_sphere_setup_tests();

  return g_test_run();
}
