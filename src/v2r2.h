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
  /* Data */
  size_t ndims;
  /* Methods */
  void (*dispose)(V2RObject *);
  void (*free)(V2RObject *);
  void (*copy)(V2RObject *, V2RObject *);
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

DllExport V2RObject *v2r_object_new(V2RObjectType *, double *, double *, double *);

#endif
