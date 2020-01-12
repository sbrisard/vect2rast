/**
 * @file v2r_spheroid.h
 *
 * @brief This header defines n-dimensional spheres.
 */
#ifndef __V2R_SPHEROID_H_202001120654__
#define __V2R_SPHEROID_H_202001120654__

#include "vect2rast.h"

#define V2R_SPHEROID_DATA(spheroid) ((V2R_SpheroidData *)((spheroid)->data))

typedef struct V2R_SpheroidData_ {
  double equatorial_radius;
  double polar_radius;
  double *axis;
  double q1, q2;
} V2R_SpheroidData;

DllExport double v2r_spheroid_equatorial_radius(V2R_Object const *spheroid);
DllExport double v2r_spheroid_polar_radius(V2R_Object const *spheroid);
DllExport void v2r_spheroid_axis(V2R_Object const *spheroid, double *axis);
DllExport V2R_Object *v2r_spheroid_new(double const *center,
                                       double equatorial_radius,
                                       double polar_radius, double const *axis);

#endif
