#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <vect2rast.h>

typedef struct {
  V2R_Object *object;
  double *point;
  bool expected;
} V2R_ObjectBelongsTestData;

V2R_ObjectBelongsTestData *v2r_object_belongs_test_data_new(V2R_Object *object,
                                                            double *point,
                                                            bool expected) {
  V2R_ObjectBelongsTestData *data = malloc(sizeof(V2R_ObjectBelongsTestData));
  data->object = v2r_object_copy(object);
  data->point = malloc(object->type->ndims * sizeof(double));
  data->expected = expected;
  for (size_t i = 0; i < object->type->ndims; i++) {
    data->point[i] = point[i];
  }
  return data;
}

void v2r_object_belongs_test_data_free(V2R_ObjectBelongsTestData *data) {
  v2r_object_free(data->object);
  free(data->point);
  free(data);
}

void test_v2r_object_belongs(V2R_ObjectBelongsTestData *data) {
  g_assert_cmpuint(data->object->type->belongs(data->object, data->point), ==,
                   data->expected);
}

void test_disk_new() {
  size_t ndims = 2;
  double c[] = {1.2, -3.4};
  double r = 7.8;

  V2R_Object *disk = v2r_disk_new(c, r);
  g_assert_cmpint(disk->type->ndims, ==, ndims);
  g_assert_cmpfloat(v2r_ndsphere_radius(disk), ==, r);

  for (size_t i = 0; i < ndims; i++) {
    g_assert_cmpfloat(disk->center[i], ==, c[i]);
    g_assert_cmpfloat(disk->bbmin[i], ==, c[i] - r);
    g_assert_cmpfloat(disk->bbmax[i], ==, c[i] + r);
  }

  v2r_object_free(disk);
}

void test_disk_setup_tests() {
  g_test_add_func("/disk/new", test_disk_new);

  char *test_name = malloc(256);

  double c[] = {1.2, -3.4};
  double r = 5.6;
  V2R_Object *disk = v2r_disk_new(c, r);

  size_t num_points = 10;
  const double r1 = 0.95 * r;
  const double r2 = 1.05 * r;

  for (size_t i = 0; i < num_points; i++) {
    double const theta = 2 * M_PI * i / (double)num_points;
    double const sin_theta = sin(theta);
    double const cos_theta = cos(theta);
    double const p1[] = {c[0] + r1 * cos_theta, c[1] + r1 * sin_theta};

    sprintf(test_name, "/disk/belongs/true/%d*pi:%d", (int)(2 * i),
            (int)num_points);
    g_test_add_data_func_full(
        test_name, v2r_object_belongs_test_data_new(disk, p1, true),
        test_v2r_object_belongs, v2r_object_belongs_test_data_free);

    double const p2[] = {c[0] + r2 * cos_theta, c[1] + r2 * sin_theta};

    sprintf(test_name, "/disk/belongs/false/%d*pi:%d", (int)(2 * i),
            (int)num_points);
    g_test_add_data_func_full(
        test_name, v2r_object_belongs_test_data_new(disk, p2, false),
        test_v2r_object_belongs, v2r_object_belongs_test_data_free);
  }

  v2r_object_free(disk);
  free(test_name);
}

void test_sphere_new() {
  size_t ndims = 3;
  double c[] = {1.2, -3.4, 5.6};
  double r = 7.8;

  V2R_Object *sphere = v2r_sphere_new(c, r);
  g_assert_cmpint(sphere->type->ndims, ==, ndims);
  g_assert_cmpfloat(v2r_ndsphere_radius(sphere), ==, r);

  for (size_t i = 0; i < ndims; i++) {
    g_assert_cmpfloat(sphere->center[i], ==, c[i]);
    g_assert_cmpfloat(sphere->bbmin[i], ==, c[i] - r);
    g_assert_cmpfloat(sphere->bbmax[i], ==, c[i] + r);
  }

  v2r_object_free(sphere);
}

void test_sphere_setup_tests() {
  g_test_add_func("/sphere/new", test_sphere_new);
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  test_disk_setup_tests();
  test_sphere_setup_tests();

  return g_test_run();
}
