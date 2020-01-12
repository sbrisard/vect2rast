#include "v2r_ndsphere.h"

V2R_NDSphereData *v2r_ndsphere_data_new(double radius) {
  V2R_NDSphereData *data = malloc(sizeof(V2R_NDSphereData));
  data->radius = radius;
  data->radius2 = radius * radius;
  return data;
}

void *v2r_ndsphere_data_copy(void const *data_) {
  V2R_NDSphereData const *data = data_;
  return v2r_ndsphere_data_new(data->radius);
}

void v2r_ndsphere_get_bounding_box(V2R_Object *sphere,
                                   double *bbmin,
                                   double *bbmax) {
  double const r = V2R_NDSPHERE_DATA(sphere)->radius;
  double *end = sphere->center + sphere->type->dim;
  for (double *x = sphere->center, *x1 = bbmin, *x2 = bbmax;
       x < end; x += 1, x1 += 1, x2 += 1) {
    *x1 = *x - r;
    *x2 = *x + r;
  }
}

double v2r_ndsphere_radius(V2R_Object const *sphere) {
  return V2R_NDSPHERE_DATA(sphere)->radius;
}

bool v2r_disk_belongs(V2R_Object const *disk, double const *point) {
  const double x = point[0] - disk->center[0];
  const double y = point[1] - disk->center[1];
  V2R_DiskData *data = disk->data;
  return x * x + y * y <= data->radius2;
}

V2R_ObjectType const Disk = {.name = "Disk",
    .dim = 2,
    .belongs = v2r_disk_belongs,
    .data_copy = v2r_ndsphere_data_copy,
    .data_free = free,
    .get_bounding_box = v2r_ndsphere_get_bounding_box};

V2R_Object *v2r_disk_new(double const *center, double radius) {
  V2R_Object *object = v2r_object_new(&Disk);

  const double x = center[0];
  const double y = center[1];

  object->center[0] = x;
  object->center[1] = y;

  object->data = v2r_ndsphere_data_new(radius);

  return object;
}

bool v2r_sphere_belongs(V2R_Object const *sphere, double const *point) {
  const double x = point[0] - sphere->center[0];
  const double y = point[1] - sphere->center[1];
  const double z = point[2] - sphere->center[2];
  V2R_SphereData *data = sphere->data;
  return x * x + y * y + z * z <= data->radius2;
}

V2R_ObjectType const Sphere = {.name = "Sphere",
    .dim = 3,
    .belongs = v2r_sphere_belongs,
    .data_copy = v2r_ndsphere_data_copy,
    .data_free = free,
    .get_bounding_box = v2r_ndsphere_get_bounding_box};

V2R_Object *v2r_sphere_new(double const *center, double radius) {
  V2R_Object *object = v2r_object_new(&Sphere);

  const double x = center[0];
  const double y = center[1];
  const double z = center[2];

  object->center[0] = x;
  object->center[1] = y;
  object->center[2] = z;
  object->data = v2r_ndsphere_data_new(radius);

  return object;
}
