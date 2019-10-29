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

typedef struct V2RSphereData_ {
  double radius;
} V2RSphereData;

#endif
