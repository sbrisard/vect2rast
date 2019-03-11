#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "voxelize.h"

typedef struct {
  Sphere *sphere;
  double *point;
  bool expected;
} SphereBelongsTestData;

SphereBelongsTestData *sphere_belongs_test_data_new(Sphere *sphere,
                                                    double *point,
                                                    bool expected) {
  SphereBelongsTestData *data = g_new(SphereBelongsTestData, 1);
  data->sphere = sphere_copy(sphere);
  data->point = g_new(double, sphere->ndims);
  data->expected = expected;
  for (size_t i = 0; i < sphere->ndims; i++) {
    data->point[i] = point[i];
  }
  return data;
}

void sphere_belongs_test_data_free(SphereBelongsTestData *data) {
  g_free(data->sphere);
  g_free(data->point);
  g_free(data);
}

void test_sphere_belongs(SphereBelongsTestData *test_data) {
  g_assert_cmpuint(
      test_data->sphere->belongs(test_data->sphere, test_data->point), ==,
      test_data->expected);
}

void test_sphere_bbox(Sphere *sphere) {
  double *min = g_new(double, sphere->ndims);
  double *max = g_new(double, sphere->ndims);
  sphere->bbox(sphere, min, max);
  for (size_t i = 0; i < sphere->ndims; i++) {
    g_assert_cmpfloat(min[i], ==, sphere->center[i] - sphere->radius);
    g_assert_cmpfloat(max[i], ==, sphere->center[i] + sphere->radius);
  }
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  double c[] = {1.2, -3.4, 5.6};
  double r = 7.8;
  Sphere *sphere = sphere_new(3, c, r);
  double theta = 0.35;
  double phi = 1.9;
  double n[] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
  double s = 0.95;
  double p1[] = {c[0] + s * r * n[0], c[1] + s * r * n[1], c[2] + s * r * n[2]};
  s = 1.05;
  double p2[] = {c[0] + s * r * n[0], c[1] + s * r * n[1], c[2] + s * r * n[2]};
  g_test_add_data_func_full("/sphere/belongs/true",
                            sphere_belongs_test_data_new(sphere, p1, true),
                            test_sphere_belongs, sphere_belongs_test_data_free);
  g_test_add_data_func_full("/sphere/belongs/false",
                            sphere_belongs_test_data_new(sphere, p2, false),
                            test_sphere_belongs, sphere_belongs_test_data_free);
  g_test_add_data_func_full("/sphere/bbox", sphere_copy(sphere),
                            test_sphere_bbox, sphere_free);
  sphere_free(sphere);
  return g_test_run();
}
