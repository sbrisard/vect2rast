/**
 * @file
 */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

#include "voxelize.h"

static void init_bounds(double x_min, double x_max, double h_inv, int *i_min,
                        int *i_max) {
  *i_min = ceil(h_inv * x_min - 0.5);
  *i_max = floor(h_inv * x_max - 0.5);
}

static void particle_free(Particle *particle) {
  g_free(particle->center);
  g_free(particle);
}

static Particle *particle_copy(Particle *, Particle *);

static void particle3d_voxelize(Particle *particle, double *dim, size_t *size,
                                guint8 *grid, guint8 value) {
  const size_t ndims = 3;
  const double L0 = dim[0], L1 = dim[1], L2 = dim[2];
  const int n0 = size[0], n1 = size[1], n2 = size[2];
  const double h0 = L0 / n0, h1 = L1 / n1, h2 = L2 / n2;

  double *min = g_new(double, ndims);
  double *max = g_new(double, ndims);
  particle->bbox(particle, min, max);
  int i0_min, i0_max, i1_min, i1_max, i2_min, i2_max;
  init_bounds(min[0], max[0], n0 / L0, &i0_min, &i0_max);
  init_bounds(min[1], max[1], n1 / L1, &i1_min, &i1_max);
  init_bounds(min[2], max[2], n2 / L2, &i2_min, &i2_max);
  g_free(max);
  g_free(min);

  double *x = g_new(double, ndims);
  int j0 = i0_min >= 0 ? i0_min : i0_min + n0;
  x[0] = (i0_min + 0.5) * h0;
  for (int i0 = i0_min; i0 <= i0_max; i0++) {
    int j1 = i1_min >= 0 ? i1_min : i1_min + n1;
    x[1] = (i1_min + 0.5) * h1;
    for (int i1 = i1_min; i1 <= i1_max; i1++) {
      int j2 = i2_min >= 0 ? i2_min : i2_min + n2;
      x[2] = (i2_min + 0.5) * h2;
      for (int i2 = i2_min; i2 <= i2_max; i2++) {
        if (particle->belongs(particle, x)) {
          grid[(j0 * n1 + j1) * n2 + j2] = value;
        }
        x[2] += h2;
        ++j2;
        if (j2 == n2) {
          j2 = 0;
        }
      }
      ++j1;
      x[1] += h1;
      if (j1 == n1) {
        j1 = 0;
      }
    }
    ++j0;
    x[0] += h0;
    if (j0 == n0) {
      j0 = 0;
    }
  }
  g_free(x);
}

static void particle2d_voxelize(Particle *particle, double *dim, size_t *size,
                                guint8 *grid, guint8 value) {
  const size_t ndims = 2;
  const double L0 = dim[0], L1 = dim[1];
  const int n0 = size[0], n1 = size[1];
  const double h0 = L0 / n0, h1 = L1 / n1;

  double *min = g_new(double, ndims);
  double *max = g_new(double, ndims);
  particle->bbox(particle, min, max);
  int i0_min, i0_max, i1_min, i1_max;
  init_bounds(min[0], max[0], n0 / L0, &i0_min, &i0_max);
  init_bounds(min[1], max[1], n1 / L1, &i1_min, &i1_max);
  g_free(max);
  g_free(min);

  double *x = g_new(double, ndims);
  int j0 = i0_min >= 0 ? i0_min : i0_min + n0;
  x[0] = (i0_min + 0.5) * h0;
  for (int i0 = i0_min; i0 <= i0_max; i0++) {
    int j1 = i1_min >= 0 ? i1_min : i1_min + n1;
    x[1] = (i1_min + 0.5) * h1;
    for (int i1 = i1_min; i1 <= i1_max; i1++) {
      if (particle->belongs(particle, x)) {
        grid[j0 * n1 + j1] = value;
      }
      ++j1;
      x[1] += h1;
      if (j1 == n1) {
        j1 = 0;
      }
    }
    ++j0;
    x[0] += h0;
    if (j0 == n0) {
      j0 = 0;
    }
  }
  g_free(x);
}

void particle_init(Particle *particle, size_t ndims, double *center,
                   particle_free_t free, particle_copy_t copy,
                   particle_belongs_t belongs, particle_bbox_t bbox) {
  particle->ndims = ndims;
  if (particle->center == NULL) {
    particle->center = g_new(double, ndims);
  }
  for (size_t i = 0; i < ndims; i++) {
    particle->center[i] = center[i];
  }
  particle->free = free;
  particle->copy = copy;
  particle->belongs = belongs;
  particle->bbox = bbox;
  if (ndims == 2) {
    particle->voxelize = particle2d_voxelize;
  } else if (ndims == 3) {
    particle->voxelize = particle3d_voxelize;
  } else {
    particle->voxelize = NULL;
  }
}

Particle *particle_new(size_t ndims, double *center) {
  Particle *particle = g_new0(Particle, 1);
  particle->free = particle_free;
  particle->copy = particle_copy;
  particle_init(particle, ndims, center, particle_free, particle_copy, NULL,
                NULL);
  return particle;
}

Particle *particle_copy(Particle *src, Particle *dest) {
  if (dest == NULL) {
    /* Use g_new0 so as to make sure that dest->center is indeed NULL, and will
     * be allocated by particle_init. */
    dest = g_new0(Particle, 1);
  }
  particle_init(dest, src->ndims, src->center, src->free, src->copy,
                src->belongs, src->bbox);
  /* TODO This does not allow for alternate implementations of the
     voxelize() method. */
  return dest;
}

static bool sphere2d_belongs(Particle *particle, double *point) {
  Sphere *sphere = (Sphere *)particle;
  double x = point[0] - sphere->center[0];
  double y = point[1] - sphere->center[1];
  double r = sphere->radius;
  return x * x + y * y <= r * r;
}

static bool sphere3d_belongs(Particle *particle, double *point) {
  Sphere *sphere = (Sphere *)particle;
  double x = point[0] - sphere->center[0];
  double y = point[1] - sphere->center[1];
  double z = point[2] - sphere->center[2];
  double r = sphere->radius;
  return x * x + y * y + z * z <= r * r;
}

static void sphere2d_bbox(Particle *particle, double *min, double *max) {
  Sphere *sphere = (Sphere *)particle;
  double r = sphere->radius;
  double *c = sphere->center;
  *min = (*c) - r;
  *max = (*c) + r;
  ++min;
  ++max;
  ++c;
  *min = (*c) - r;
  *max = (*c) + r;
}

static void sphere3d_bbox(Particle *particle, double *min, double *max) {
  Sphere *sphere = (Sphere *)particle;
  double r = sphere->radius;
  double *c = sphere->center;
  *min = (*c) - r;
  *max = (*c) + r;
  ++min;
  ++max;
  ++c;
  *min = (*c) - r;
  *max = (*c) + r;
  ++min;
  ++max;
  ++c;
  *min = (*c) - r;
  *max = (*c) + r;
}

static Particle *sphere_copy(Particle *, Particle *);

Sphere *sphere_new(size_t ndims, double *center, double radius) {
  Sphere *sphere = g_new0(Sphere, 1);
  particle_belongs_t sphere_belongs;
  particle_bbox_t sphere_bbox;
  if (ndims == 2) {
    sphere_belongs = sphere2d_belongs;
    sphere_bbox = sphere2d_bbox;
  } else if (ndims == 3) {
    sphere_belongs = sphere3d_belongs;
    sphere_bbox = sphere3d_bbox;
  } else {
    sphere_belongs = NULL;
    sphere_bbox = NULL;
  }
  particle_init(sphere, ndims, center, particle_free, sphere_copy,
                sphere_belongs, sphere_bbox);
  sphere->radius = radius;
  return sphere;
}

Particle *sphere_copy(Particle *src, Particle *dest) {
  Sphere *dest_ = dest == NULL ? g_new0(Sphere, 1) : dest;
  particle_copy(src, dest_);
  dest_->radius = ((Sphere *)src)->radius;
  return dest_;
}

static void spheroid_bbox(Particle *particle, double *min, double *max) {
  Spheroid *spheroid = (Spheroid *)particle;
  const double a = spheroid->equatorial_radius;
  const double c = spheroid->polar_radius;
  const double a2 = a * a;
  const double c2_m_a2 = c * c - a2;
  for (size_t i = 0; i < spheroid->ndims; i++) {
    const double d_i = spheroid->axis[i];
    const double dx_i = sqrt(a2 + c2_m_a2 * d_i * d_i);
    min[i] = spheroid->center[i] - dx_i;
    max[i] = spheroid->center[i] + dx_i;
  }
}

static bool spheroid_belongs(Particle *particle, double *point) {
  Spheroid *spheroid = (Spheroid *)particle;
  double x_dot_x = 0.;
  double d_dot_x = 0.;
  for (size_t i = 0; i < particle->ndims; i++) {
    const double x_i = point[i] - spheroid->center[i];
    x_dot_x += x_i * x_i;
    d_dot_x += spheroid->axis[i] * x_i;
  }
  return spheroid->_q1 * x_dot_x + spheroid->_q2 * d_dot_x * d_dot_x <= 1.;
}

particle_copy_t spheroid_copy;
/* static Particle *spheroid_copy(Particle *, Particle *); */

/*
Spheroid *spheroid_new(size_t ndims, double *center, double equatorial_radius,
                       double polar_radius, double *axis) {
  Spheroid *spheroid = g_new(Spheroid, 1);
  spheroid->ndims = ndims;
  spheroid->center = g_new(double, ndims);
  spheroid->equatorial_radius = equatorial_radius;
  spheroid->polar_radius = polar_radius;
  spheroid->axis = g_new(double, ndims);
  for (size_t i = 0; i < ndims; i++) {
    spheroid->center[i] = center[i];
    spheroid->axis[i] = axis[i];
  }
  spheroid->_q1 = 1. / (equatorial_radius * equatorial_radius);
  spheroid->_q2 = 1. / (polar_radius * polar_radius) - spheroid->_q1;
  spheroid->copy = spheroid_copy;
  spheroid->belongs = spheroid_belongs;
  spheroid->bbox = spheroid_bbox;
  spheroid->voxelize = particle3d_voxelize;
  return spheroid;
}

void spheroid_free(Spheroid *spheroid) {
  g_free(spheroid->center);
  g_free(spheroid->axis);
  g_free(spheroid);
}

Particle *spheroid_copy(Particle *particle) {
  Spheroid *spheroid = (Spheroid *)particle;
  return spheroid_new(spheroid->ndims, spheroid->center,
                      spheroid->equatorial_radius, spheroid->polar_radius,
                      spheroid->axis);
}
*/
