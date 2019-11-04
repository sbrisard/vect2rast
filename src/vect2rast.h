#ifndef __V2R_H__
#define __V2R_H__

#include <stdbool.h>
#include <stdlib.h>

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

typedef struct V2R_Object_ V2R_Object;
typedef struct V2R_ObjectType_ V2R_ObjectType;

struct V2R_ObjectType_ {
  char *name;
  size_t ndims;
  size_t data_size;

  bool (*belongs)(V2R_Object const *, double const *);
};

struct V2R_Object_ {
  V2R_ObjectType *type;
  double *center;
  double *bbmin;
  double *bbmax;
  void *data;
};

DllExport V2R_Object *v2r_object_new(V2R_ObjectType const *);
DllExport void v2r_object_free(V2R_Object *);
DllExport V2R_Object *v2r_object_copy(V2R_Object const *);

DllExport double v2r_ndsphere_radius(V2R_Object const *);
DllExport V2R_Object *v2r_disk_new(double const *, double);
DllExport V2R_Object *v2r_sphere_new(double const *, double);

#endif
