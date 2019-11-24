#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vect2rast.h"

V2R_Object *v2r_object_new(V2R_ObjectType const *type) {
  size_t const dim = type->dim;
  V2R_Object *object = malloc(sizeof(V2R_Object));
  object->type = type;

  size_t const size = dim * sizeof(double);
  object->center = malloc(size);
  object->bbmin = malloc(size);
  object->bbmax = malloc(size);

  object->data = NULL;

  return object;
}

void v2r_object_free(V2R_Object *object) {
  object->type->data_free(object->data);
  free(object->center);
  free(object->bbmin);
  free(object->bbmax);
}

V2R_Object *v2r_object_copy(V2R_Object *const object) {
  V2R_Object *copy = v2r_object_new(object->type);
  size_t const dim = object->type->dim;
  for (size_t i = 0; i < dim; i++) {
    copy->center[i] = object->center[i];
    copy->bbmin[i] = object->bbmin[i];
    copy->bbmax[i] = object->bbmax[i];
  }
  copy->data = object->type->data_copy(object->data);
  return copy;
}
