#ifndef FLUIDINITIALIZE_H
#define FLUIDINITIALIZE_H

#include "global_definitions.h"
#include "FluidQuantites.h"
#include "Vector3D.h"

class FluidInitialize
{
public:
  FluidInitialize() = default;
  void operator() (const Vector3DR &x, const Vector3DR &y, const Vector3DR &z,
                   FluidQuantity3DArray *q) noexcept;
};

#endif // FLUIDINITIALIZE_H
