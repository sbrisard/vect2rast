/**
 * @file vect2rast.h
 *
 * @brief This header defines the basic object model, as well as the
 * raterization functions.
 */

#ifndef __VECT2RAST_H_202001120825__
#define __VECT2RAST_H_202001120825__

#include <stdbool.h>
#include <stdlib.h>

#include <string>

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

typedef struct V2R_Object_ V2R_Object;

typedef struct V2R_ObjectType_ V2R_ObjectType;

/**
 * @brief The base structure that defines the type of a geometrical object.
 *
 * Defines some "class variables" and the object's "methods".
 *
 * @see V2R_ObjectType
 */
struct V2R_ObjectType_ {
  /** @brief A `NULL`-terminated string that holds the name of the type. */
  std::string name;
  /** @brief The dimension of the embedding space (usually 2 or 3). */
  size_t dim;

  /**
   * @brief Return a copy of a geometric object's (type-specific) `data` field.
   *
   * This function is required by ::v2r_object_copy.
   */
  void *(*data_copy)(void const *);

  /**
   * @brief Deallocate a geometric object's (type-specific) `data` field.
   *
   * This function is required by ::v2r_object_free.
   */
  void (*data_free)(void *);

  /** @brief Return true if the specified `point` belongs to `object`. */
  bool (*belongs)(V2R_Object const *object, double const *point);

  /**
   * @brief Update `bbmin` and `bbmax` with a bounding box to `object`.
   *
   * The bounding box must satisfy the condition that if, for some `point`,
   * there exists an index `i` in the `0, ..., object->type->dim` range
   * such that `point[i] < bbmin[i]` or `point[i] > bbmax[i]`, then
   * `object->type->belongs(object, point)` must return `false`.
   *
   * The bounding box is not required to be optimal. However, tighter bounding
   * boxes should result in faster rasterization.
   */
  void (*get_bounding_box)(V2R_Object const *object, double *bbmin, double *bbmax);
};

/**
 * @brief The base structure that defines a geometrical object.
 *
 * @see V2R_Object
 */
struct V2R_Object_ {
  /** @brief The type of the geometrical object (the container that holds its
   * methods). */
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

/**
 * @brief Allocate a new geometric object of the specified `type`.
 *
 * This function also allocates the V2R_Object_::center array.
 */
DllExport V2R_Object *v2r_object_new(V2R_ObjectType *type);

/**
 * @brief Create a copy of the specified geometric `object`.
 *
 * The returned copy is located at the specified `center`. If `center == NULL`,
 * then the copy is centered at the same location as the original objec.
 *
 * This function requires V2R_ObjectType_::data_copy.
 */
DllExport V2R_Object *v2r_object_copy(V2R_Object const *object,
                                      double const *center);

/**
 * @brief Deallocate the specified geometric `object`.
 *
 * This function requires V2R_ObjectType_::data_free.
 */
DllExport void v2r_object_free(V2R_Object *object);

DllExport int v2r_raster(V2R_Object const *object, double const *length,
                         size_t const *size, int *grid, int value);

#endif
