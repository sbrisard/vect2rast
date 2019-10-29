#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "v2r2.h"

V2RObject *v2r_object_new(V2RObjectType *type, double *center, double *bbmin,
                          double *bbmax) {
  const size_t ndims = type->ndims;
  V2RObject *object = malloc(sizeof(V2RObject));

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


int main(int argc, char** argv) {
  printf("coucou\n");
  return 0;
}
