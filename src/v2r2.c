#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

#include "v2r2.h"

V2RObject *v2r_object_new(V2RObjectType *type, double *center, double *bbmin,
                          double *bbmax) {
  const size_t ndims = type->ndims;
  V2RObject *object = (V2RObject *)g_new(V2RObject, 1);

  const size_t size = ndims * sizeof(double);
  object->center = malloc(size);
  object->bbmin = malloc(size);
  object->bbmax = malloc(size);

  object->data = NULL;

  for (size_t i = 0; i < object->type->ndims; i++) {
    object->center[i] = center[i];
    object->bbmin[i] = bbmin[i];
    object->bbmax[i] = bbmax[i];
  }

  return object;
}

void v2r_object_free(V2RObject *object) {
  if (object->type->dispose) {
    object->type->dispose(object);
  }
  free(object->center);
  free(object->bbmin);
  free(object->bbmax);
}
