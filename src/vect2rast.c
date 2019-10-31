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

struct SphereData_ {
  double radius;
  double radius2;
};

typedef struct SphereData_ DiskData;
typedef struct SphereData_ SphereData;

SphereData *v2r_sphere_data_new(double radius) {
  SphereData *data = malloc(sizeof(SphereData));
  data->radius = radius;
  data->radius2 = radius * radius;
  return data;
}


#define V2R_NDSPHERE_DATA(sphere) ((SphereData *)((sphere)->data))

void v2r_sphere_data_free(SphereData *data) { free(data); }

bool v2r_disk_belongs(V2RObject *disk, double *point) {
  const double x = point[0] - disk->center[0];
  const double y = point[1] - disk->center[1];
  DiskData *data = disk->data;
  return x * x + y * y <= data->radius2;
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
  object->bbmin[0] = x - radius;
  object->bbmax[0] = x + radius;

  object->center[1] = y;
  object->bbmin[1] = y - radius;
  object->bbmax[1] = y + radius;

  object->data = v2r_sphere_data_new(radius);

  return object;
}

double v2r_ndsphere_radius(V2RObject const *sphere){
  return V2R_NDSPHERE_DATA(sphere)->radius;
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

  printf("r = %g\n", V2R_NDSPHERE_DATA(disk)->radius);
  printf("r**2 = %g\n", V2R_NDSPHERE_DATA(disk)->radius2);

  return 0;
}
