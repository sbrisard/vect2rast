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

typedef struct Particle Particle;
typedef void particle_free_t(Particle *);
typedef Particle *particle_copy_t(Particle *, Particle *);
typedef bool particle_belongs_t(Particle *, double *);
typedef void particle_voxelize_t(Particle *, double *, size_t *, guint8 *,
                                    guint8);

struct Particle {
  /* Data */
  size_t ndims;
  double *center;
  double *bbmin;
  double *bbmax;
  /* Methods */
  particle_free_t *free;
  particle_copy_t *copy;
  particle_belongs_t *belongs;
  particle_voxelize_t *voxelize;
};

DllExport Particle *particle_new(size_t, double *, double *, double *);

typedef struct Sphere {
  struct Particle;
  double radius;
} Sphere;

DllExport Sphere *sphere_new(size_t ndims, double *center, double radius);

typedef struct Spheroid {
  struct Particle;
  double equatorial_radius;
  double polar_radius;
  double *axis;
  double _q1, _q2;
} Spheroid;

DllExport Spheroid *spheroid_new(double *, double, double, double *);

#endif
