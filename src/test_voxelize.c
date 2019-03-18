#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "voxelize.h"

#define minimum_image(x, L)  \
  {                          \
    double L_half = .5 * L;  \
    if (x < -L_half) x += L; \
    if (x > L_half) x -= L;  \
  }

void init_icosahedron(double *vertex) {
  const double phi = .5 * (1. + sqrt(5.));
  const double u = 1. / sqrt(1 + phi * phi);
  const double v = phi * u;
  vertex[0] = 0.;
  vertex[1] = -u;
  vertex[2] = -v;
  vertex += 3;
  vertex[0] = 0.;
  vertex[1] = -u;
  vertex[2] = +v;
  vertex += 3;
  vertex[0] = 0.;
  vertex[1] = +u;
  vertex[2] = -v;
  vertex += 3;
  vertex[0] = 0.;
  vertex[1] = +u;
  vertex[2] = +v;
  vertex += 3;
  vertex[0] = -u;
  vertex[1] = -v;
  vertex[2] = 0.;
  vertex += 3;
  vertex[0] = -u;
  vertex[1] = +v;
  vertex[2] = 0.;
  vertex += 3;
  vertex[0] = +u;
  vertex[1] = -v;
  vertex[2] = 0.;
  vertex += 3;
  vertex[0] = +u;
  vertex[1] = +v;
  vertex[2] = 0.;
  vertex += 3;
  vertex[0] = -v;
  vertex[1] = 0.;
  vertex[2] = -u;
  vertex += 3;
  vertex[0] = -v;
  vertex[1] = 0.;
  vertex[2] = +u;
  vertex += 3;
  vertex[0] = +v;
  vertex[1] = 0.;
  vertex[2] = -u;
  vertex += 3;
  vertex[0] = +v;
  vertex[1] = 0.;
  vertex[2] = +u;
  vertex += 3;
}

typedef struct {
  Particle *particle;
  double *point;
  bool expected;
} ParticleBelongsTestData;

ParticleBelongsTestData *particle_belongs_test_data_new(Particle *particle,
                                                        double *point,
                                                        bool expected) {
  ParticleBelongsTestData *data = g_new(ParticleBelongsTestData, 1);
  data->particle = particle->copy(particle, NULL);
  data->point = g_new(double, particle->ndims);
  data->expected = expected;
  for (size_t i = 0; i < particle->ndims; i++) {
    data->point[i] = point[i];
  }
  return data;
}

void particle_belongs_test_data_free(ParticleBelongsTestData *data) {
  g_free(data->particle);
  g_free(data->point);
  g_free(data);
}

void test_sphere_belongs(ParticleBelongsTestData *test_data) {
  Particle *particle = test_data->particle;
  g_assert_cmpuint(particle->belongs(particle, test_data->point), ==,
                   test_data->expected);
}

void test_sphere_bbox(Sphere *sphere) {
  for (size_t i = 0; i < sphere->ndims; i++) {
    g_assert_cmpfloat(sphere->bbmin[i], ==, sphere->center[i] - sphere->radius);
    g_assert_cmpfloat(sphere->bbmax[i], ==, sphere->center[i] + sphere->radius);
  }
}

typedef struct {
  Particle *particle;
  double *dim;
  size_t *size;
} ParticleVoxelizeTestData;

ParticleVoxelizeTestData *particle_voxelize_test_data_new(Particle *particle,
                                                          double *dim,
                                                          size_t *size) {
  size_t ndims = particle->ndims;
  ParticleVoxelizeTestData *data = g_new(ParticleVoxelizeTestData, 1);
  data->particle = particle->copy(particle, NULL);
  data->dim = g_new(double, ndims);
  data->size = g_new(size_t, ndims);
  for (size_t i = 0; i < ndims; i++) {
    data->dim[i] = dim[i];
    data->size[i] = size[i];
  }
  return data;
}

void particle_voxelize_test_data_free(ParticleVoxelizeTestData *data) {
  /* TODO: this is not a polymorphic call. */
  data->particle->free(data->particle);
  g_free(data->dim);
  g_free(data->size);
}

void test_particle_voxelize(ParticleVoxelizeTestData *data) {
  const size_t ndims = data->particle->ndims;

  size_t n[] = {data->size[0], data->size[1], 1};
  double L[] = {data->dim[0], data->dim[1], 0.};
  double c[] = {data->particle->center[0], data->particle->center[1], 0.};
  if (ndims == 3) {
    const size_t i = ndims - 1;
    n[i] = data->size[i];
    L[i] = data->dim[i];
    c[i] = data->particle->center[i];
  }
  double h[3];
  for (size_t i = 0; i < 3; i++) {
    h[i] = L[i] / n[i];
  }

  guint8 *actual = g_new0(guint8, n[0] * n[1] * n[2]);

  data->particle->voxelize(data->particle, data->dim, data->size, actual, 1);

  /* Create copy of particle, centered at the origin, since we will compute the
   * minimum image of the CP vector (C: center; P: current point). */
  Particle *particle = data->particle->copy(data->particle, NULL);
  for (size_t i = 0; i < ndims; i++) {
    particle->center[i] = 0.;
  }
  double *point = g_new(double, ndims);
  for (size_t i0 = 0; i0 < n[0]; i0++) {
    point[0] = (i0 + 0.5) * h[0] - c[0];
    minimum_image(point[0], L[0]);
    for (size_t i1 = 0; i1 < n[1]; i1++) {
      point[1] = (i1 + 0.5) * h[1] - c[1];
      minimum_image(point[1], L[1]);
      for (size_t i2 = 0; i2 < n[2]; i2++) {
        point[2] = (i2 + 0.5) * h[2] - c[2];
        minimum_image(point[2], L[2]);
        size_t j = (i0 * n[1] + i1) * n[2] + i2;
        guint8 expected = particle->belongs(particle, point);
        g_assert_cmpuint(expected, ==, actual[j]);
      }
    }
  }
  particle->free(particle);
  g_free(point);
  g_free(actual);
}

void setup_sphere_tests() {
  size_t ndims = 3;
  double c[] = {1.2, -3.4, 5.6};
  double r = 7.8;
  Sphere *sphere = sphere_new(3, c, r);
  char *name = g_new(char *, 255);

  /* Sphere->bbmin, Sphere->bbmax */
  g_test_add_data_func_full("/sphere/bbox", sphere->copy(sphere, NULL),
                            test_sphere_bbox, sphere->free);

  /* Sphere->belongs() */
  double s1 = 0.95;
  double s2 = 1.05;
  size_t ndirs = 12;
  double *n = g_new(double, ndims *ndirs);
  init_icosahedron(n);
  double *p = g_new(double, ndims);
  for (size_t j = 0; j < ndirs; j++) {
    for (size_t i = 0; i < ndims; i++) {
      p[i] = c[i] + s1 * r * n[ndims * j + i];
    }
    sprintf(name, "/sphere/belongs/true/%d", j);
    g_test_add_data_func_full(
        name, particle_belongs_test_data_new(sphere, p, true),
        test_sphere_belongs, particle_belongs_test_data_free);
    sprintf(name, "/sphere/belongs/false/%d", j);
    for (size_t i = 0; i < ndims; i++) {
      p[i] = c[i] + s2 * r * n[ndims * j + i];
    }
    g_test_add_data_func_full(
        name, particle_belongs_test_data_new(sphere, p, false),
        test_sphere_belongs, particle_belongs_test_data_free);
  }
  g_free(p);
  g_free(n);
  sphere->free(sphere);

  /* Sphere->voxelize() */
  double dim[] = {1.5, 2.6, 3.7};
  size_t size[] = {50, 60, 70};
  double xi = 0.05;
  double eta = 0.06;
  double zeta = 0.07;
  r = 0.5;

  for (size_t i0 = 0; i0 <= 1; i0++) {
    c[0] = i0 == 0 ? xi * dim[0] : (1. - xi) * dim[0];
    for (size_t i1 = 0; i1 <= 1; i1++) {
      c[1] = i1 == 0 ? eta * dim[1] : (1. - eta) * dim[1];
      for (size_t i2 = 0; i2 <= 1; i2++) {
        c[2] = i2 == 0 ? zeta * dim[2] : (1. - zeta) * dim[2];
        sprintf(name, "/sphere/voxelize/%d", (i0 * 2 + i1) * 2 + i2);
        sphere = sphere_new(ndims, c, r);
        g_test_add_data_func_full(
            name, particle_voxelize_test_data_new(sphere, dim, size),
            test_particle_voxelize, particle_voxelize_test_data_free);
        sphere->free(sphere);
      }
    }
  }

  g_free(name);
}

void setup_spheroid_tests() {
  double c[] = {0.75, 1.3, 1.85};
  double dim[] = {1.5, 2.6, 3.7};
  size_t size[] = {50, 60, 70};

  double theta = 0.3;
  double phi = 0.5;
  double d[] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};

  Spheroid *spheroid = spheroid_new(c, 0.37, .19, d);
  ParticleVoxelizeTestData *data =
      particle_voxelize_test_data_new(spheroid, dim, size);
  spheroid->free(spheroid);

  g_test_add_data_func_full("/spheroid/voxelize/1", data,
                            test_particle_voxelize,
                            particle_voxelize_test_data_free);
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  setup_sphere_tests();
  setup_spheroid_tests();

  return g_test_run();
}
