#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "voxelize.h"

#define minimum_image(x, L)  \
  {                          \
    double L_half = .5 * L;  \
    if (x < -L_half) x += L; \
    if (x > L_half) x -= L;  \
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
  /* TODO: memory leak. */
  particle->free(particle);
  g_free(point);
  g_free(actual);
}

void setup_sphere_tests() {
  double c[] = {1.2, -3.4, 5.6};
  double r = 7.8;
  Sphere *sphere = sphere_new(3, c, r);
  double theta = 0.35;
  double phi = 1.9;
  double n[] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
  double s = 0.95;
  double p1[] = {c[0] + s * r * n[0], c[1] + s * r * n[1], c[2] + s * r * n[2]};
  s = 1.05;
  double p2[] = {c[0] + s * r * n[0], c[1] + s * r * n[1], c[2] + s * r * n[2]};

  ParticleBelongsTestData *data1, *data2;
  data1 = particle_belongs_test_data_new(sphere, p1, true);
  data2 = particle_belongs_test_data_new(sphere, p2, false);
  g_test_add_data_func_full("/sphere/belongs/1", data1, test_sphere_belongs,
                            particle_belongs_test_data_free);
  g_test_add_data_func_full("/sphere/belongs/2", data2, test_sphere_belongs,
                            particle_belongs_test_data_free);
  Sphere *sphere2 = sphere->copy(sphere, NULL);
  g_test_add_data_func_full("/sphere/bbox", sphere2, test_sphere_bbox,
                            sphere2->free);
  sphere->free(sphere);

  c[0] = 0.75;
  c[1] = 1.3;
  c[2] = 1.85;
  r = 0.69;
  double dim[] = {1.5, 2.6, 3.7};
  size_t size[] = {50, 60, 70};

  sphere = sphere_new(3, c, r);
  ParticleVoxelizeTestData *data3;
  data3 = particle_voxelize_test_data_new(sphere, dim, size);
  sphere->free(sphere);

  g_test_add_data_func_full("/sphere/voxelize/1", data3, test_particle_voxelize,
                            particle_voxelize_test_data_free);
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
