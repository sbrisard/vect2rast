#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vect2rast.h>

V2RObject *v2r_object_new(V2RObjectType const *type) {
  const size_t ndims = type->ndims;
  V2RObject *object = malloc(sizeof(V2RObject));
  object->type = type;

  const size_t size = ndims * sizeof(double);
  object->center = malloc(size);
  object->bbmin = malloc(size);
  object->bbmax = malloc(size);

  object->data = NULL;

  return object;
}

void v2r_object_free(V2RObject *object) {
  free(object->center);
  free(object->bbmin);
  free(object->bbmax);
  free(object->data);
}

V2RObject *v2r_object_copy(V2RObject const *object) {
  V2RObject *copy = v2r_object_new(object->type);
  for (size_t i = 0; i < object->type->ndims; i++) {
    copy->center[i] = object->center[i];
    copy->bbmin[i] = object->bbmin[i];
    copy->bbmax[i] = object->bbmax[i];
  }
  copy->data = malloc(object->type->data_size);
  memcpy(copy->data, object->data, object->type->data_size);
  return copy;
}

/* V2RObjectType const Sphere3D = {.name = "Sphere", */
/*                                 .ndims = 3, */
/*                                 .data_size = 2 * sizeof(double), */
/*                                 .belongs = v2r_sphere_belongs}; */
