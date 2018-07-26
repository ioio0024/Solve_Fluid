#ifndef REFLECTBOUNDARY_H
#define REFLECTBOUNDARY_H

#include "global_definitions.h"
#include "FluidQuantites.h"

struct ReflectBoundaryXLeft
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct ReflectBoundaryXRight
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct ReflectBoundaryYLeft
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct ReflectBoundaryYRight
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct ReflectBoundaryZLeft
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};

struct ReflectBoundaryZRight
{
  void operator()(Int32 is, Int32 ie,
                  Int32 js, Int32 je,
                  Int32 ks, Int32 ke,
                  FluidQuantity3DArray *sia,
                  FluidQuantity3DArray *via
                 ) noexcept;
};
#endif // REFLECTBOUNDARY_H
