#ifndef IS_BOUNDARY_CONDITION_H
#define IS_BOUNDARY_CONDITION_H

#include <utility>
#include "global_definitions.h"

template <typename T>
class Vector3D;

struct is_boundary_condtion_impl
{
  template <typename F>
  static auto check(F*) -> decltype(
      std::declval<F>(Int32 /*is*/, Int32 /*ie*/,
                      Int32 /*js*/, Int32 /*je*/,
                      Int32 /*ks*/, Int32 /*ke*/,
                      Vector3D<Real> */*d_sia*/,
                      Vector3D<Real> */*u_sia*/, Vector3D<Real> */*v_sia*/, Vector3D<Real> */*w_sia*/,
                      Vector3D<Real> */*p_sia*/,
                      Vector3D<Real> */*d_via*/,
                      Vector3D<Real> */*u_via*/, Vector3D<Real> */*v_via*/, Vector3D<Real> */*w_via*/,
                      Vector3D<Real> */*p_via*/),
      std::true_type());
};

#endif // IS_BOUNDARY_CONDITION_H
