/**
 * @file
 */

#ifndef __VOXELIZE_H__
#define __VOXELIZE_H__

#include <stdlib.h>

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

typedef struct Sphere_ {
  /** Number of spatial dimensions **/
  size_t ndims;
  /** Coordinates of center **/
  double *center;
  /** Radius **/
  double radius;
  /** Returns `true` if point belongs to sphere. **/
  bool (*belongs)(struct Sphere_ *, double *);
} Sphere;

DllExport Sphere *sphere_new(size_t ndims, double *center, double radius);
DllExport void sphere_free(Sphere *);
DllExport Sphere *sphere_copy(Sphere *);
DllExport bool sphere_belongs(Sphere *, double *);
DllExport void sphere_bbox(Sphere *, double *, double *);

#endif
