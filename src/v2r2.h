#ifndef __V2R_H__
#define __V2R_H__

#include <stdlib.h>

#if _WIN32
#define DllExport __declspec(dllexport)
#else
#define DllExport
#endif

typedef struct V2RObject_ V2RObject;
typedef struct V2RObjectType_ V2RObjectType;

struct V2RObjectType_ {
  size_t ndims;
  size_t data_size;

  void (*free)(V2RObject *);
  V2RObject * (*copy)(V2RObject *);
  bool (*belongs)(V2RObject *, double *);
};

struct V2RObject_ {
  V2RObjectType *type;
  double *center;
  double *bbmin;
  double *bbmax;
  void *data;
};

typedef struct V2RSphereData_ {
  double radius;
} V2RSphereData;

DllExport V2RObject *v2r_object_new(V2RObjectType *);

#endif
