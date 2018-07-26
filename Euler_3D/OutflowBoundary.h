#ifndef OUTFLOWBOUNDARY_H
#define OUTFLOWBOUNDARY_H

#include "global_definitions.h"
#include "FluidQuantites.h"

struct OutflowBoundaryXLeft
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct OutflowBoundaryXRight
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct OutflowBoundaryYLeft
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct OutflowBoundaryYRight
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct OutflowBoundaryZLeft
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct OutflowBoundaryZRight
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};
#endif // OUTFLOWBOUNDARY_H
