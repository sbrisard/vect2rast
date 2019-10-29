#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "v2r2.h"

V2RObject *v2r_object_new(V2RObjectType const *type) {
  const size_t ndims = type->ndims;
  V2RObject *object = malloc(sizeof(V2RObject));

  const size_t size = ndims * sizeof(double);
  object->center = malloc(size);
  object->bbmin = malloc(size);
  object->bbmax = malloc(size);

  object->data = malloc(type->data_size);

  return object;
}

void v2r_object_free(V2RObject *object) {
  free(object->center);
  free(object->bbmin);
  free(object->bbmax);
  free(object->data);
}

V2RObject *v2r_object_copy(V2RObject *object) {
  V2RObject *copy = v2r_object_new(object->type);
  for (size_t i = 0; i < object->type->ndims; i++) {
    copy->center[i] = object->center[i];
    copy->bbmin[i] = object->bbmin[i];
    copy->bbmax[i] = object->bbmax[i];
  }
  memcpy(copy->data, object->data, object->type->data_size);
  return copy;
}

#define V2R_SPHERE_RADIUS_INDEX 0
#define V2R_SPHERE_SQR_RADIUS_INDEX 1

bool v2r_sphere_belongs(V2RObject *sphere, double *point) {
  double r2 = 0.0;
  for (size_t i = 0; i < sphere->type->ndims; i++) {
    const double x_i = point[i] - sphere->center[i];
    r2 += x_i * x_i;
  }
  return r2 <= V2R_OBJECT_DOUBLE_AT(sphere, V2R_SPHERE_SQR_RADIUS_INDEX);
}

V2RObjectType const Sphere2D = {.ndims = 2,
                                .data_size = 2 * sizeof(double),
                                .belongs = v2r_sphere_belongs};

V2RObjectType const Sphere3D = {.ndims = 3,
                                .data_size = 2 * sizeof(double),
                                .belongs = v2r_sphere_belongs};

V2RObject *v2r_sphere_new(size_t ndims, double *center, double radius) {
  V2RObject *object;
  if (ndims == 2) {
    object = v2r_object_new(&Sphere2D);
  } else if (ndims == 3) {
    object = v2r_object_new(&Sphere3D);
  } else {
    return NULL;
  }

  for (size_t i = 0; i < ndims; i++) {
    const double x_i = center[i];
    object->center[i] = x_i;
    object->bbmin[i] = x_i - radius;
    object->bbmax[i] = x_i + radius;
  }

  V2R_OBJECT_DOUBLE_AT(object, V2R_SPHERE_RADIUS_INDEX) = radius;
  V2R_OBJECT_DOUBLE_AT(object, V2R_SPHERE_SQR_RADIUS_INDEX) = radius * radius;

  return object;
}

int main(int argc, char **argv) {
  printf("coucou\n");

  double center[] = {1., 2.};
  double radius = 0.5;
  V2RObject *sphere = v2r_sphere_new(2, center, radius);

  printf("r**2 = %g\n",
         V2R_OBJECT_DOUBLE_AT(sphere, V2R_SPHERE_SQR_RADIUS_INDEX));

  return 0;
}
