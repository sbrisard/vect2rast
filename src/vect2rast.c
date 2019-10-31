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

  object->data = malloc(type->data_size);

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
  memcpy(copy->data, object->data, object->type->data_size);
  return copy;
}

bool v2r_disk_belongs(V2RObject *disk, double *point) {
  double r2 = 0.0;
  const double x = point[0] - disk->center[0];
  const double y = point[1] - disk->center[1];
  return x * x + y * y <=
         V2R_OBJECT_DOUBLE_AT(disk, V2R_SPHERE_SQR_RADIUS_INDEX);
}

V2RObjectType const Disk = {.name = "Disk",
                            .ndims = 2,
                            .data_size = 2 * sizeof(double),
                            .belongs = v2r_disk_belongs};

V2RObject *v2r_disk_new(double *center, double radius) {
  V2RObject *object = v2r_object_new(&Disk);

  const double x = center[0];
  const double y = center[1];

  object->center[0] = x;
  object->bbmin[0] = x-radius;
  object->bbmax[0] = x+radius;

  object->center[1] = y;
  object->bbmin[1] = y-radius;
  object->bbmax[1] = y+radius;

  V2R_OBJECT_DOUBLE_AT(object, V2R_SPHERE_RADIUS_INDEX) = radius;
  V2R_OBJECT_DOUBLE_AT(object, V2R_SPHERE_SQR_RADIUS_INDEX) = radius * radius;

  return object;
}

/* V2RObjectType const Sphere3D = {.name = "Sphere", */
/*                                 .ndims = 3, */
/*                                 .data_size = 2 * sizeof(double), */
/*                                 .belongs = v2r_sphere_belongs}; */


int main(int argc, char **argv) {
  printf("coucou\n");

  double center[] = {1., 2.};
  double radius = 0.5;
  V2RObject *disk = v2r_disk_new(center, radius);

  printf("Type name = %s\n", disk->type->name);

  printf("r**2 = %g\n",
         V2R_OBJECT_DOUBLE_AT(disk, V2R_SPHERE_SQR_RADIUS_INDEX));

  return 0;
}
