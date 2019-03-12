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

typedef struct Particle {
  size_t ndims;
  double *center;
  struct Particle *(*copy)(struct Particle *);
  bool (*belongs)(struct Particle *, double *);
  void (*bbox)(struct Particle *, double *, double *);
  void (*voxelize)(struct Particle *, double *, size_t *, guint8 *, guint8);
} Particle;

DllExport Particle *particle_new(size_t ndims, double *center);
DllExport void particle_free(Particle *);

typedef struct {
  struct Particle;
  double radius;
} Sphere;

DllExport Sphere *sphere_new(size_t ndims, double *center, double radius);
DllExport void sphere_free(Sphere *);

typedef struct {
  struct Particle;
  double equatorial_radius;
  double polar_radius;
  double *axis;
  double _q1, _q2;
} Spheroid;

DllExport Spheroid *spheroid_new(size_t, double *, double, double, double *);
DllExport void spheroid_free(Spheroid *);

#endif
