#include <vect2rast.h>

#define V2R_NDSPHERE_DATA(sphere) ((SphereData *)((sphere)->data))

struct SphereData_ {
  double radius;
  double radius2;
};

typedef struct SphereData_ DiskData;
typedef struct SphereData_ SphereData;

void v2r_sphere_data_free(SphereData *data) { free(data); }

SphereData *v2r_sphere_data_new(double radius) {
  SphereData *data = malloc(sizeof(SphereData));
  data->radius = radius;
  data->radius2 = radius * radius;
  return data;
}

bool v2r_disk_belongs(V2RObject *disk, double *point) {
  const double x = point[0] - disk->center[0];
  const double y = point[1] - disk->center[1];
  DiskData *data = disk->data;
  return x * x + y * y <= data->radius2;
}

V2RObjectType const Disk = {.name = "Disk",
                            .ndims = 2,
                            .data_size = sizeof(DiskData),
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
