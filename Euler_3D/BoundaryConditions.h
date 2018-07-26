#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
::BoundaryConditions() noexcept
{
  BoundaryConditionsX();
  BoundaryConditionsY();
  BoundaryConditionsZ();
}

#include "BoundaryConditionsX.h"
#include "BoundaryConditionsY.h"
#include "BoundaryConditionsZ.h"

#endif // BOUNDARYCONDITIONS_H
