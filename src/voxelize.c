#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

#include "voxelize.h"

static particle_copy_t particle_copy;
static particle_copy_t sphere_copy;
static particle_copy_t spheroid_copy;

static void init_bounds(double x_min, double x_max, double h_inv, int *i_min,
                        int *i_max) {
  *i_min = ceil(h_inv * x_min - 0.5);
  *i_max = floor(h_inv * x_max - 0.5);
}

static void *particle_allocate(Particle *particle, size_t ndims) {
  particle->ndims = ndims;
  particle->center = g_new(double, ndims);
  particle->bbmin = g_new(double, ndims);
  particle->bbmax = g_new(double, ndims);
  return particle;
}

static void particle_free(Particle *particle) {
  g_free(particle->center);
  g_free(particle);
}

static void particle3d_voxelize(Particle *particle, double *dim, size_t *size,
                                guint8 *grid, guint8 value) {
  const size_t ndims = 3;
  const double L0 = dim[0], L1 = dim[1], L2 = dim[2];
  const int n0 = size[0], n1 = size[1], n2 = size[2];
  const double h0 = L0 / n0, h1 = L1 / n1, h2 = L2 / n2;

  int i0_min, i0_max, i1_min, i1_max, i2_min, i2_max;
  init_bounds(particle->bbmin[0], particle->bbmax[0], n0 / L0, &i0_min,
              &i0_max);
  init_bounds(particle->bbmin[1], particle->bbmax[1], n1 / L1, &i1_min,
              &i1_max);
  init_bounds(particle->bbmin[2], particle->bbmax[2], n2 / L2, &i2_min,
              &i2_max);

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

  int i0_min, i0_max, i1_min, i1_max;
  init_bounds(particle->bbmin[0], particle->bbmax[0], n0 / L0, &i0_min,
              &i0_max);
  init_bounds(particle->bbmin[1], particle->bbmax[1], n1 / L1, &i1_min,
              &i1_max);

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

Particle *particle_new(size_t ndims, double *center, double *bbmin,
                       double *bbmax) {
  Particle *particle = g_new(Particle, 1);
  particle_allocate(particle, ndims);
  for (size_t i = 0; i < ndims; i++) {
    particle->center[i] = center[i];
    particle->bbmin[i] = bbmin[i];
    particle->bbmax[i] = bbmax[i];
  }
  particle->free = particle_free;
  particle->copy = particle_copy;
  if (ndims == 2) {
    particle->voxelize = particle2d_voxelize;
  } else if (ndims == 3) {
    particle->voxelize = particle3d_voxelize;
  } else {
    particle->voxelize = NULL;
  }
  return particle;
}

Particle *particle_copy(Particle *src, Particle *dest) {
  if (dest == NULL) {
    dest = g_new(Particle, 1);
    particle_allocate(dest, src->ndims);
  }
  for (size_t i = 0; i < src->ndims; i++) {
    dest->center[i] = src->center[i];
    dest->bbmin[i] = src->bbmin[i];
    dest->bbmax[i] = src->bbmax[i];
  }
  dest->free = src->free;
  dest->copy = src->copy;
  dest->belongs = src->belongs;
  dest->voxelize = src->voxelize;
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

Sphere *sphere_new(size_t ndims, double *center, double radius) {
  Sphere *sphere = g_new(Sphere, 1);
  particle_allocate(sphere, ndims);
  for (size_t i = 0; i < ndims; i++) {
    sphere->center[i] = center[i];
    sphere->bbmin[i] = center[i] - radius;
    sphere->bbmax[i] = center[i] + radius;
  }
  sphere->radius = radius;
  sphere->free = particle_free;
  sphere->copy = sphere_copy;
  if (ndims == 2) {
    sphere->belongs = sphere2d_belongs;
    sphere->voxelize = particle2d_voxelize;
  } else if (ndims == 3) {
    sphere->belongs = sphere3d_belongs;
    sphere->voxelize = particle3d_voxelize;
  } else {
    sphere->belongs = NULL;
    sphere->voxelize = NULL;
  }
  return sphere;
}

Particle *sphere_copy(Particle *src, Particle *dest) {
  Sphere *src_ = src;
  Sphere *dest_ = dest;
  if (dest_ == NULL) {
    dest_ = g_new(Sphere, 1);
    particle_allocate(dest_, src->ndims);
  }
  particle_copy(src, dest_);
  dest_->radius = src_->radius;
  return dest_;
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

void spheroid_allocate(Spheroid *spheroid) {
  particle_allocate(spheroid, 3);
  spheroid->axis = g_new(double, spheroid->ndims);
}

void spheroid_free(Spheroid *spheroid) {
  g_free(spheroid->axis);
  particle_free(spheroid);
}

Spheroid *spheroid_new(double *center, double equatorial_radius,
                       double polar_radius, double *axis) {
  Spheroid *spheroid = g_new(Spheroid, 1);
  spheroid->ndims = 3;
  spheroid_allocate(spheroid);
  spheroid->equatorial_radius = equatorial_radius;
  spheroid->polar_radius = polar_radius;
  spheroid->_q1 = 1. / (equatorial_radius * equatorial_radius);
  spheroid->_q2 = 1. / (polar_radius * polar_radius) - spheroid->_q1;
  const double a2 = equatorial_radius * equatorial_radius;
  const double c2_m_a2 = polar_radius * polar_radius - a2;
  for (size_t i = 0; i < spheroid->ndims; i++) {
    /* TODO Normalize axis? */
    spheroid->center[i] = center[i];
    spheroid->axis[i] = axis[i];
    const double r = sqrt(a2 + c2_m_a2 * axis[i] * axis[i]);
    spheroid->bbmin[i] = spheroid->center[i] - r;
    spheroid->bbmax[i] = spheroid->center[i] + r;
  }
  spheroid->free = spheroid_free;
  spheroid->copy = spheroid_copy;
  spheroid->belongs = spheroid_belongs;
  spheroid->voxelize = particle3d_voxelize;
  return spheroid;
}

Particle *spheroid_copy(Particle *src, Particle *dest) {
  Spheroid *src_ = src;
  Spheroid *dest_ = dest;
  if (dest_ == NULL) {
    dest_ = g_new(Spheroid, 1);
    spheroid_allocate(dest_);
  }
  particle_copy(src, dest_);
  dest_->equatorial_radius = src_->equatorial_radius;
  dest_->polar_radius = src_->polar_radius;
  dest_->_q1 = src_->_q1;
  dest_->_q2 = src_->_q2;
  for (size_t i = 0; i < src_->ndims; i++) {
    dest_->axis[i] = src_->axis[i];
  }
  return dest_;
}
