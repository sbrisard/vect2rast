#include <glib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "v2r.h"

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

void test_particle_belongs(ParticleBelongsTestData *test_data) {
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
        test_particle_belongs, particle_belongs_test_data_free);
    sprintf(name, "/sphere/belongs/false/%d", j);
    for (size_t i = 0; i < ndims; i++) {
      p[i] = c[i] + s2 * r * n[ndims * j + i];
    }
    g_test_add_data_func_full(
        name, particle_belongs_test_data_new(sphere, p, false),
        test_particle_belongs, particle_belongs_test_data_free);
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

void fill_basis(double *e1, double *e2, double *e3) {
  const double x = e1[0];
  const double y = e1[1];
  const double z = e1[2];
  const double s = 1. / hypot(y, z);
  e2[0] = 0.;
  e2[1] = s * z;
  e2[2] = -s * y;
  e3[0] = -1. / s;
  e3[1] = s * x * y;
  e3[2] = s * x * z;
}

void test_spheroid_bbox(Sphere *spheroid) {
  /* TODO: implement test. */
}

void setup_spheroid_bbox_tests() {
  char *name = g_new(char, 255);
  const size_t ndims = 3;
  const size_t ndirs = 12;
  double c[] = {1.2, -3.4, 5.6};

  double *n = g_new(double, ndims *ndirs);
  init_icosahedron(n);

  Spheroid *spheroid;
  size_t test_index = 0;
  for (size_t i = 0; i < ndirs; i++) {
    spheroid = spheroid_new(c, 0.78, 0.19, n + i * ndims);
    spheroid->free(spheroid);
    sprintf(name, "/spheroid/oblate/bbox/%d", test_index);

    spheroid = spheroid_new(c, 0.19, 0.78, n + i * ndims);
    spheroid->free(spheroid);
    sprintf(name, "/spheroid/prolate/bbox/%d", test_index);

    ++test_index;
  }

  g_free(name);
  g_free(n);
}

void setup_spheroid_belongs_tests() {
  char *name = g_new(char, 255);
  const size_t ndims = 3;
  const size_t ndirs = 12;
  double c[] = {1.2, -3.4, 5.6};

  double *n = g_new(double, ndims *ndirs);
  init_icosahedron(n);
  double *e1 = g_new(double, ndims);
  double *e2 = g_new(double, ndims);
  double *p = g_new(double, ndims);
  double *radii = g_new(double, ndims);

  double s1 = 0.95;
  double s2 = 1.05;
  double r1 = 0.78;
  double r2 = 0.19;
  double x1, x2, x3, x, y, z;

  for (size_t i = 0; i < ndirs; i++) {
    double *e3 = n + i * ndims;
    fill_basis(e3, e1, e2);
    Spheroid *oblate = spheroid_new(c, 0.78, 0.19, e3);
    Spheroid *prolate = spheroid_new(c, 0.19, 0.78, e3);
    for (size_t j = 0; j < ndirs; j++) {
      double *n_j = n + ndims * j;
      x1 = n_j[0] * oblate->equatorial_radius;
      x2 = n_j[1] * oblate->equatorial_radius;
      x3 = n_j[2] * oblate->polar_radius;
      x = x1 * e1[0] + x2 * e2[0] + x3 * e3[0];
      y = x1 * e1[1] + x2 * e2[1] + x3 * e3[1];
      z = x1 * e1[2] + x2 * e2[2] + x3 * e3[2];

      p[0] = c[0] + s1 * x;
      p[1] = c[1] + s1 * y;
      p[2] = c[2] + s1 * z;
      sprintf(name, "/spheroid/oblate/belongs/true/%d", i * ndirs + j);
      g_test_add_data_func_full(
          name, particle_belongs_test_data_new(oblate, p, true),
          test_particle_belongs, particle_belongs_test_data_free);
      p[0] = c[0] + s2 * x;
      p[1] = c[1] + s2 * y;
      p[2] = c[2] + s2 * z;
      sprintf(name, "/spheroid/oblate/belongs/false/%d", i * ndirs + j);
      g_test_add_data_func_full(
          name, particle_belongs_test_data_new(oblate, p, false),
          test_particle_belongs, particle_belongs_test_data_free);
    }
    oblate->free(oblate);
    prolate->free(prolate);
  }

  g_free(name);
  g_free(n);
  g_free(e1);
  g_free(e2);
  g_free(p);
}

void setup_spheroid_voxelize_tests() {
  size_t ndims = 3;
  size_t ndirs = 12;
  double *n = g_new(double, ndims *ndirs);
  init_icosahedron(n);
  double dim[] = {1.5, 2.6, 3.7};
  size_t size[] = {50, 60, 70};
  double xi = 0.05;
  double eta = 0.06;
  double zeta = 0.07;

  double *c;
  c = g_new(double, ndims);
  Spheroid *spheroid;
  size_t test_index = 0;
  char *name = g_new(char, 255);
  for (size_t i0 = 0; i0 <= 1; i0++) {
    c[0] = i0 == 0 ? xi * dim[0] : (1. - xi) * dim[0];
    for (size_t i1 = 0; i1 <= 1; i1++) {
      c[1] = i1 == 0 ? eta * dim[1] : (1. - eta) * dim[1];
      for (size_t i2 = 0; i2 <= 1; i2++) {
        c[2] = i2 == 0 ? zeta * dim[2] : (1. - zeta) * dim[2];
        for (size_t j = 0; j < ndirs; j++) {
          spheroid = spheroid_new(c, 0.78, 0.19, n + j * ndims);
          sprintf(name, "/spheroid/oblate/voxelize/%d", test_index);
          g_test_add_data_func_full(
              name, particle_voxelize_test_data_new(spheroid, dim, size),
              test_particle_voxelize, particle_voxelize_test_data_free);
          spheroid->free(spheroid);
          spheroid = spheroid_new(c, 0.19, 0.78, n + j * ndims);
          sprintf(name, "/spheroid/prolate/voxelize/%d", test_index);
          g_test_add_data_func_full(
              name, particle_voxelize_test_data_new(spheroid, dim, size),
              test_particle_voxelize, particle_voxelize_test_data_free);
          spheroid->free(spheroid);
          ++test_index;
        }
      }
    }
  }
  g_free(n);
  g_free(c);
  g_free(name);
}

void setup_spheroid_tests() {
  setup_spheroid_bbox_tests();
  setup_spheroid_belongs_tests();
  setup_spheroid_voxelize_tests();
}

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  setup_sphere_tests();
  setup_spheroid_tests();

  return g_test_run();
}
