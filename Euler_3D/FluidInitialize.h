#ifndef FLUIDINITIALIZE_H
#define FLUIDINITIALIZE_H

#include "global_definitions.h"
#include "Vector3D.h"

class FluidInitialize
{
public:
  FluidInitialize() = default;
  void operator() (const Rvec &x, const Rvec &y, const Rvec &z,
                   Vector3DR *d, Vector3DR *u, Vector3DR *v, Vector3DR *w,
                   Vector3DR *p) noexcept;
};

#endif // FLUIDINITIALIZE_H
