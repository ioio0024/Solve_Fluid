#ifndef FLUIDINITIALIZE_H
#define FLUIDINITIALIZE_H

#include "global_definitions.h"

void FluidInitialize(const Rvec& x,
                     Rvec *density, Rvec *velocity, Rvec *pressure);
#endif // FLUIDINITIALIZE_H
