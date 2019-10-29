#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "v2r2.h"

V2RObject *v2r_object_new(V2RObjectType *type) {
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

double v2r_sphere_radius2(V2RObject *sphere) {
  return ((double *)sphere->data)[1];
}

bool v2r_sphere2d_belongs(V2RObject *sphere, double *point) {
  double x = point[0] - sphere->center[0];
  double y = point[1] - sphere->center[1];
  return x * x + y * y <= v2r_sphere_radius2(sphere);
}

const V2RObjectType Sphere2D = {.ndims = 2,
				.free = v2r_object_free,
				.copy = NULL,
				.belongs = v2r_sphere2d_belongs};

V2RObject *v2r_sphere2d_new(double *center, double radius) {
  V2RObject *object = v2r_object_new(&Sphere2D);

  const double x0 = center[0];
  object->center[0] = x0;
  object->bbmin[0] = x0 - radius;
  object->bbmax[0] = x0 + radius;

  const double x1 = center[1];
  object->center[1] = x1;
  object->bbmin[1] = x1 - radius;
  object->bbmax[1] = x1 + radius;

  double *data = malloc(2 * sizeof(double));
  data[0] = radius;
  data[1] = radius * radius;
  object->data = data;

  return object;
}


int main(int argc, char **argv) {
  printf("coucou\n");


  double center[] = {1., 2.};
  double radius = 0.5;
  V2RObject *sphere = v2r_sphere2d_new(center, radius);

  printf("r**2 = %g\n", v2r_sphere_radius2(sphere));

  return 0;
}
