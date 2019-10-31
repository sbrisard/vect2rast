#ifndef __V2R_H__
#define __V2R_H__

#include <stdlib.h>

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

#define V2R_OBJECT_DOUBLE_AT(object, i) ((double *)(object)->data)[i]

typedef struct V2RObject_ V2RObject;
typedef struct V2RObjectType_ V2RObjectType;

struct V2RObjectType_ {
  char *name;
  size_t ndims;
  size_t data_size;

  bool (*belongs)(V2RObject *, double *);
};

struct V2RObject_ {
  V2RObjectType *type;
  double *center;
  double *bbmin;
  double *bbmax;
  void *data;
};

DllExport V2RObject *v2r_object_new(V2RObjectType const *);
DllExport void v2r_object_free(V2RObject *);
DllExport V2RObject *v2r_object_copy(V2RObject const *);

#define V2R_SPHERE_RADIUS_INDEX 0
#define V2R_SPHERE_SQR_RADIUS_INDEX 1

DllExport V2RObject *v2r_disk_new(double *, double);

#endif
