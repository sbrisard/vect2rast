/**
 * @file v2r_ndsphere.h
 *
 * @brief This header defines n-dimensional spheres.
 */
#ifndef __V2R_NDSPHERE_H_202001120646__
#define __V2R_NDSPHERE_H_202001120646__

#include "vect2rast.h"

#define V2R_NDSPHERE_DATA(sphere) ((V2R_NDSphereData *)((sphere)->data))

struct V2R_NDSphereData_ {
  double radius;
  double radius2;
};

typedef struct V2R_NDSphereData_ V2R_NDSphereData;
typedef struct V2R_NDSphereData_ V2R_DiskData;
typedef struct V2R_NDSphereData_ V2R_SphereData;

DllExport double v2r_ndsphere_radius(V2R_Object const *sphere);
DllExport V2R_Object *v2r_disk_new(double const *center, double radius);
DllExport V2R_Object *v2r_sphere_new(double const *center, double radius);

#endif
