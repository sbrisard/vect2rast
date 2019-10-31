#include <vect2rast.h>

#define V2R_NDSPHERE_DATA(sphere) ((V2R_NDSphereData *)((sphere)->data))

struct V2R_NDSphereData_ {
  double radius;
  double radius2;
};

typedef struct V2R_NDSphereData_ V2R_NDSphereData;
typedef struct V2R_NDSphereData_ V2R_DiskData;
typedef struct V2R_NDSphereData_ V2R_SphereData;

void v2r_ndsphere_data_free(V2R_NDSphereData *data) { free(data); }

V2R_NDSphereData *v2r_ndsphere_data_new(double radius) {
  V2R_NDSphereData *data = malloc(sizeof(V2R_NDSphereData));
  data->radius = radius;
  data->radius2 = radius * radius;
  return data;
}

bool v2r_disk_belongs(V2R_Object *disk, double *point) {
  const double x = point[0] - disk->center[0];
  const double y = point[1] - disk->center[1];
  V2R_DiskData *data = disk->data;
  return x * x + y * y <= data->radius2;
}

V2R_ObjectType const Disk = {.name = "Disk",
                            .ndims = 2,
                            .data_size = sizeof(V2R_DiskData),
                            .belongs = v2r_disk_belongs};

V2R_Object *v2r_disk_new(double *center, double radius) {
  V2R_Object *object = v2r_object_new(&Disk);

  const double x = center[0];
  const double y = center[1];

  object->center[0] = x;
  object->bbmin[0] = x - radius;
  object->bbmax[0] = x + radius;

  object->center[1] = y;
  object->bbmin[1] = y - radius;
  object->bbmax[1] = y + radius;

  object->data = v2r_ndsphere_data_new(radius);

  return object;
}

double v2r_ndsphere_radius(V2R_Object const *sphere){
  return V2R_NDSPHERE_DATA(sphere)->radius;
}
