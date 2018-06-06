#ifndef FLUIDINITIALIZE_H
#define FLUIDINITIALIZE_H

#include "global_definitions.h"
#include "Vector2D.h"

class FluidInitialize
{
public:
  FluidInitialize() = default;
  void operator() (const Rvec &x, const Rvec &y,
                   Vector2D *d, Vector2D *u, Vector2D *v, Vector2D *p) noexcept;
};

#endif // FLUIDINITIALIZE_H
