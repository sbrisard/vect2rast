/**
 * @file vect2rast.h
 *
 * Vector-to-raster conversion of simple geometrical objects
 */

#ifndef __V2R_H__
#define __V2R_H__

#include <stdbool.h>
#include <stdlib.h>

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

/**
 * @see V2R_Object_
 */
typedef struct V2R_Object_ V2R_Object;

/**
 * @see V2R_ObjectType_
 */
typedef struct V2R_ObjectType_ V2R_ObjectType;

/**
 * @brief The base structure that defines the type of a geometrical object.
 *
 * Defines some "class variables" and the object's "methods".
 */
struct V2R_ObjectType_ {
  /** @brief A `NULL`-terminated string that holds the name of the type. */
  char *name;
  /** @brief The dimension of the embedding space (usually 2 or 3). */
  size_t dim;

  void *(*data_copy)(void const *);
  void (*data_free)(void *);
  bool (*belongs)(V2R_Object const *object, double const *point);
  void (*get_bounding_box)(V2R_Object *object, double *bbmin, double *bbmax);
};

/**
 * @brief The base structure that defines a geometrical object.
 */
struct V2R_Object_ {
  /** @brief The type of the geometrical object. */
  V2R_ObjectType *type;
  /**
   * @brief The Coordinates of the center.
   *
   * Array of length `object->type->dim`.
   *
   * @see V2R_ObjectType_::dim.
   */
  double *center;
  /**
   * @brief Additional data that may define the shape, orientation, ...
   *
   * The actual structure of the data depends on the type of the object.
   */
  void *data;
};

DllExport V2R_Object *v2r_object_new(V2R_ObjectType const *object);
DllExport V2R_Object *v2r_object_copy(V2R_Object const *object,
                                      double const *center);
DllExport void v2r_object_free(V2R_Object *object);

DllExport int v2r_raster(V2R_Object *object, double *length, size_t *size,
                         int *grid, int value);

DllExport double v2r_ndsphere_radius(V2R_Object const *sphere);
DllExport V2R_Object *v2r_disk_new(double const *center, double radius);
DllExport V2R_Object *v2r_sphere_new(double const *center, double radius);

DllExport double v2r_spheroid_equatorial_radius(V2R_Object const *spheroid);
DllExport double v2r_spheroid_polar_radius(V2R_Object const *spheroid);
DllExport void v2r_spheroid_axis(V2R_Object const *spheroid, double *axis);
DllExport V2R_Object *v2r_spheroid_new(double const *center,
                                       double equatorial_radius,
                                       double polar_radius, double const *axis);

#endif
