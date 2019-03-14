#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "voxelize.h"

typedef struct {
  Particle *particle;
  double *point;
  bool expected;
} ParticleBelongsTestData;

ParticleBelongsTestData *particle_belongs_test_data_new(Particle *particle,
                                                        double *point,
                                                        bool expected) {
  ParticleBelongsTestData *data = g_new(ParticleBelongsTestData, 1);
  data->particle = particle->copy(particle);
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
  double *min = g_new(double, sphere->ndims);
  double *max = g_new(double, sphere->ndims);
  sphere->bbox(sphere, min, max);
  for (size_t i = 0; i < sphere->ndims; i++) {
    g_assert_cmpfloat(min[i], ==, sphere->center[i] - sphere->radius);
    g_assert_cmpfloat(max[i], ==, sphere->center[i] + sphere->radius);
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
  data->particle = particle->copy(particle);
  data->dim = g_new(double, ndims);
  data->size = g_new(size_t, ndims);
  for (size_t i = 0; i < ndims; i++) {
    data->dim[i] = dim[i];
    data->size[i] = size[i];
  }
  return data;
}

void particle_voxelize_test_data_free(ParticleVoxelizeTestData *test_data) {
  /* TODO: this is not a polymorphic call. */
  particle_free(test_data->particle);
  g_free(test_data->dim);
  g_free(test_data->size);
}

void test_sphere_voxelize(ParticleVoxelizeTestData *data) {
  const size_t ndims = data->particle->ndims;

  const size_t n0 = data->size[0];
  const size_t n1 = data->size[1];
  const size_t n2 = ndims == 3 ? data->size[2] : 1;

  const double L0 = data->dim[0];
  const double L1 = data->dim[1];
  const double L2 = data->dim[2];

  const double h0 = L0 / n0;
  const double h1 = L1 / n1;
  const double h2 = ndims == 3 ? L2 / n2 : 0.;

  const double c0 = data->particle->center[0];
  const double c1 = data->particle->center[1];
  const double c2 = ndims == 3 ? data->particle->center[2] : 0.;

  const double r = ((Sphere *)data->particle)->radius;
  const double sqr_r = r * r;

  guint8 *actual = g_new0(guint8, n0 * n1 * n2);

  data->particle->voxelize(data->particle, data->dim, data->size, actual, 1);

  for (size_t i0 = 0; i0 < n0; i0++) {
    double x0 = (i0 + 0.5) * h0 - c0;
    if (x0 < -.5 * L0) x0 += L0;
    if (x0 > .5 * L0) x0 -= L0;
    for (size_t i1 = 0; i1 < n1; i1++) {
      double x1 = (i1 + 0.5) * h1 - c1;
      if (x1 < -.5 * L1) x1 += L1;
      if (x1 > .5 * L1) x1 -= L1;
      for (size_t i2 = 0; i2 < n2; i2++) {
        double x2 = (i2 + 0.5) * h2 - c2;
        if (x2 < -.5 * L2) x2 += L2;
        if (x2 > .5 * L2) x2 -= L2;

        size_t j = (i0 * n1 + i1) * n2 + i2;
        guint8 expected = (x0 * x0 + x1 * x1 + x2 * x2) <= sqr_r ? 1 : 0;
        g_assert_cmpuint(expected, ==, actual[j]);
      }
    }
  }
  g_free(actual);
}

void setup_sphere_voxelize_tests() {
  double c[] = {0.75, 1.3, 1.85};
  double r = 0.69;
  double dim[] = {1.5, 2.6, 3.7};
  size_t size[] = {50, 60, 70};

  Sphere *sphere = sphere_new(3, c, r);
  ParticleVoxelizeTestData *data =
      particle_voxelize_test_data_new(sphere, dim, size);
  sphere_free(sphere);

  g_test_add_data_func_full("/sphere/voxelize/1", data, test_sphere_voxelize,
                            particle_voxelize_test_data_free);
}

void test_spheroid_voxelize(ParticleVoxelizeTestData *data) {
  const size_t ndims = data->particle->ndims;
  const Spheroid *spheroid = (Spheroid *)data->particle;

  const size_t nx = data->size[0];
  const size_t ny = data->size[1];
  const size_t nz = ndims == 3 ? data->size[2] : 1;

  const double Lx = data->dim[0];
  const double Ly = data->dim[1];
  const double Lz = data->dim[2];

  const double hx = Lx / nx;
  const double hy = Ly / ny;
  const double hz = ndims == 3 ? Lz / nz : 0.;

  const double cx = data->particle->center[0];
  const double cy = data->particle->center[1];
  const double cz = ndims == 3 ? data->particle->center[2] : 0.;

  double q[ndims * ndims];
  const double q1 =
      1. / (spheroid->equatorial_radius * spheroid->equatorial_radius);
  const double q2 = 1. / (spheroid->polar_radius * spheroid->polar_radius) - q1;
  for (size_t i = 0; i < ndims; i++) {
    for (size_t j = 0; j < ndims; j++) {
      q[ndims * i + j] = q2 * spheroid->axis[i] * spheroid->axis[j];
    }
    q[ndims * i + i] += q1;
  }

  guint8 *actual = g_new0(guint8, nx * ny * nz);

  data->particle->voxelize(data->particle, data->dim, data->size, actual, 1);

  for (size_t i = 0; i < nx; i++) {
    double x = (i + 0.5) * hx - cx;
    if (x < -.5 * Lx) x += Lx;
    if (x > .5 * Lx) x -= Lx;
    for (size_t j = 0; j < ny; j++) {
      double y = (j + 0.5) * hy - cy;
      if (y < -.5 * Ly) y += Ly;
      if (y > .5 * Ly) y -= Ly;
      for (size_t k = 0; k < nz; k++) {
        double z = (k + 0.5) * hz - cz;
        if (z < -.5 * Lz) z += Lz;
        if (z > .5 * Lz) z -= Lz;

        double q_dot_cp[] = {q[0] * x + q[1] * y + q[2] * z,
                             q[3] * x + q[4] * y + q[5] * z,
                             q[6] * x + q[7] * y + q[8] * z};
        double cp_dot_q_dot_cp =
            q_dot_cp[0] * x + q_dot_cp[1] * y + q_dot_cp[2] * z;
        guint8 expected = cp_dot_q_dot_cp <= 1. ? 1 : 0;
        g_assert_cmpuint(expected, ==, actual[(i * ny + j) * nz + k]);
      }
    }
  }
  g_free(actual);
}

void setup_spheroid_voxelize_tests() {
  double c[] = {0.75, 1.3, 1.85};
  double dim[] = {1.5, 2.6, 3.7};
  size_t size[] = {50, 60, 70};

  double theta = 0.3;
  double phi = 0.5;
  double d[] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};

  Spheroid *spheroid = spheroid_new(3, c, 0.37, .19, d);
  ParticleVoxelizeTestData *data =
      particle_voxelize_test_data_new(spheroid, dim, size);
  spheroid_free(spheroid);

  g_test_add_data_func_full("/spheroid/voxelize/1", data,
                            test_spheroid_voxelize,
                            particle_voxelize_test_data_free);
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

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
  g_test_add_data_func_full(
      "/sphere/belongs/true", particle_belongs_test_data_new(sphere, p1, true),
      test_sphere_belongs, particle_belongs_test_data_free);
  g_test_add_data_func_full("/sphere/belongs/false",
                            particle_belongs_test_data_new(sphere, p2, false),
                            test_sphere_belongs,
                            particle_belongs_test_data_free);
  g_test_add_data_func_full("/sphere/bbox", sphere->copy(sphere),
                            test_sphere_bbox, sphere_free);
  sphere_free(sphere);

  setup_sphere_voxelize_tests();
  setup_spheroid_voxelize_tests();

  return g_test_run();
}
