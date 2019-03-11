/**
 * @file
 */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

#include "voxelize.h"

/** Specialization of sphere_belongs() to two-dimensional spheres (disks). */
static bool sphere2d_belongs(Particle *particle, double *point) {
  Sphere *sphere = (Sphere *)particle;
  double x = point[0] - sphere->center[0];
  double y = point[1] - sphere->center[1];
  double r = sphere->radius;
  return x * x + y * y <= r * r;
}

/** Specialization of sphere_belongs() to three-dimensional spheres (disks). */
static bool sphere3d_belongs(Particle *particle, double *point) {
  Sphere *sphere = (Sphere *)particle;
  double x = point[0] - sphere->center[0];
  double y = point[1] - sphere->center[1];
  double z = point[2] - sphere->center[2];
  double r = sphere->radius;
  return x * x + y * y + z * z <= r * r;
}

void sphere2d_bbox(Particle *particle, double *min, double *max) {
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

void sphere3d_bbox(Particle *particle, double *min, double *max) {
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

/** Create new sphere. */
Sphere *sphere_new(size_t ndims, double *center, double radius) {
  Sphere *sphere = g_new(Sphere, 1);
  sphere->ndims = ndims;
  sphere->center = g_new(double, ndims);
  for (size_t i = 0; i < ndims; i++) {
    sphere->center[i] = center[i];
  }
  sphere->radius = radius;
  if (ndims == 2) {
    sphere->belongs = sphere2d_belongs;
    sphere->bbox = sphere2d_bbox;
  } else if (ndims == 3) {
    sphere->belongs = sphere3d_belongs;
    sphere->bbox = sphere3d_bbox;
  }
  return sphere;
}

void sphere_free(Sphere *sphere) {
  g_free(sphere->center);
  g_free(sphere);
}

Sphere *sphere_copy(Sphere *sphere) {
  return sphere_new(sphere->ndims, sphere->center, sphere->radius);
}

void init_bounds(double x_min, double x_max, double h_inv, int *i_min,
                 int *i_max) {
  *i_min = ceil(h_inv * x_min - 0.5);
  *i_max = floor(h_inv * x_max - 0.5);
}

void discretize3d(double *center, double radius, double *dim, size_t *size,
                  guint8 *grid) {
  const double c0 = center[0], c1 = center[1], c2 = center[2];
  const double L0 = dim[0], L1 = dim[1], L2 = dim[2];
  const int n0 = size[0], n1 = size[1], n2 = size[2];
  const double h0 = L0 / n0, h1 = L1 / n1, h2 = L2 / n2;
  const double sqr_radius = radius * radius;

  int i0_min, i0_max, i1_min, i1_max, i2_min, i2_max;
  init_bounds(c0 - radius, c0 + radius, n0 / L0, &i0_min, &i0_max);
  init_bounds(c1 - radius, c1 + radius, n1 / L1, &i1_min, &i1_max);
  init_bounds(c2 - radius, c2 + radius, n2 / L2, &i2_min, &i2_max);

  int j0 = i0_min >= 0 ? i0_min : i0_min + n0;
  double x0 = (i0_min + 0.5) * h0 - c0;
  double sqr_x0 = x0 * x0;
  for (int i0 = i0_min; i0 <= i0_max; i0++) {
    int j1 = i1_min >= 0 ? i1_min : i1_min + n1;
    double x1 = (i1_min + 0.5) * h1 - c1;
    double sqr_x1 = x1 * x1;
    for (int i1 = i1_min; i1 <= i1_max; i1++) {
      int j2 = i2_min >= 0 ? i2_min : i2_min + n2;
      double x2 = (i2_min + 0.5) * h2 - c2;
      double sqr_x2 = x2 * x2;
      for (int i2 = i2_min; i2 <= i2_max; i2++) {
        if (sqr_x0 + sqr_x1 + sqr_x2 <= sqr_radius) {
          ++grid[(j0 * n1 + j1) * n2 + j2];
        }
        ++j2;
        x2 += h2;
        if (j2 == n2) {
          j2 = 0;
        }
        sqr_x2 = x2 * x2;
      }
      ++j1;
      x1 += h1;
      if (j1 == n1) {
        j1 = 0;
      }
      sqr_x1 = x1 * x1;
    }
    ++j0;
    x0 += h0;
    if (j0 == n0) {
      j0 = 0;
    }
    sqr_x0 = x0 * x0;
  }
}

void discretize2d(double *center, double radius, double *dim, size_t *size,
                  guint8 *grid) {
  const double c0 = center[0], c1 = center[1];
  const double L0 = dim[0], L1 = dim[1];
  const int n0 = size[0], n1 = size[1];
  const double h0 = L0 / n0, h1 = L1 / n1;
  const double sqr_radius = radius * radius;

  int i0_min, i0_max, i1_min, i1_max;
  init_bounds(c0 - radius, c0 + radius, n0 / L0, &i0_min, &i0_max);
  init_bounds(c1 - radius, c1 + radius, n1 / L1, &i1_min, &i1_max);

  int j0 = i0_min >= 0 ? i0_min : i0_min + n0;
  double x0 = (i0_min + 0.5) * h0 - c0;
  double sqr_x0 = x0 * x0;
  for (int i0 = i0_min; i0 <= i0_max; i0++) {
    int j1 = i1_min >= 0 ? i1_min : i1_min + n1;
    double x1 = (i1_min + 0.5) * h1 - c1;
    double sqr_x1 = x1 * x1;
    for (int i1 = i1_min; i1 <= i1_max; i1++) {
      if (sqr_x0 + sqr_x1 <= sqr_radius) {
        ++grid[j0 * n1 + j1];
      }
      ++j1;
      x1 += h1;
      if (j1 == n1) {
        j1 = 0;
      }
      sqr_x1 = x1 * x1;
    }
    ++j0;
    x0 += h0;
    if (j0 == n0) {
      j0 = 0;
    }
    sqr_x0 = x0 * x0;
  }
}

void test() {
  int ndims = 3;
  double dim[] = {5., 6., 7.};
  size_t size[] = {100, 200, 300};
  size_t num_spheres = 20;

  const size_t n0 = size[0], n1 = size[1], n2 = size[2];
  const double L0 = dim[0], L1 = dim[1], L2 = dim[2];
  const double h0 = L0 / n0, h1 = L1 / n1, h2 = L2 / n2;

  guint8 *expected = g_new(guint8, n0 * n1 * n2);
  guint8 *actual = g_new(guint8, n0 * n1 * n2);

  GRand *rand_ = g_rand_new_with_seed(20190224);
  double *center = g_new(double, ndims *num_spheres);
  double *radius = g_new(double, num_spheres);
  for (size_t k = 0; k < num_spheres; k++) {
    center[3 * k + 0] = g_rand_double_range(rand_, 0., L0);
    center[3 * k + 1] = g_rand_double_range(rand_, 0., L1);
    center[3 * k + 2] = g_rand_double_range(rand_, 0., L2);
    radius[k] = g_rand_double_range(rand_, 0., 0.5);
    discretize3d(center + 3 * k, radius[k], dim, size, actual);
  }

  double x0, x1, x2;
  size_t num_differences = 0;
  for (size_t i0 = 0; i0 < n0; i0++) {
    for (size_t i1 = 0; i1 < n1; i1++) {
      for (size_t i2 = 0; i2 < n2; i2++) {
        size_t j = (i0 * n1 + i1) * n2 + i2;
        expected[j] = 0;
        for (size_t k = 0; k < num_spheres; k++) {
          x0 = (i0 + 0.5) * h0 - center[3 * k + 0];
          if (x0 < -.5 * L0) x0 += L0;
          if (x0 > .5 * L0) x0 -= L0;
          x1 = (i1 + 0.5) * h1 - center[3 * k + 1];
          if (x1 < -.5 * L1) x1 += L1;
          if (x1 > .5 * L1) x1 -= L1;
          x2 = (i2 + 0.5) * h2 - center[3 * k + 2];
          if (x2 < -.5 * L2) x2 += L2;
          if (x2 > .5 * L2) x2 -= L2;
          if ((x0 * x0 + x1 * x1 + x2 * x2) <= radius[k] * radius[k]) {
            ++expected[j];
          }
        }
        g_assert_cmpuint(expected[j], ==, actual[j]);
      }
    }
  }
  g_assert_cmpuint(num_differences, ==, 0);

  g_free(center);
  g_free(radius);
  g_rand_free(rand_);
  g_free(actual);
  g_free(expected);
}
