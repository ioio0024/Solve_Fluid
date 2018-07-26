template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>::Advection()
{
  // PV to LIA x
  AdvectionX(ny_ + 1, nz_ + 1, pv_, lia_x_, &pv_new_, &lia_x_new_);
  // LIA y to SIA z
  AdvectionX(ny_, nz_ + 1, lia_y_, sia_z_, &lia_y_new_, &sia_z_new_);
  // LIA z to SIA y
  AdvectionX(ny_ + 1, nz_, lia_z_, sia_y_, &lia_z_new_, &sia_y_new_);
  // SIA x to VIA
  AdvectionX(ny_, nz_, sia_x_, via_, &sia_x_new_, &via_new_);
  UpdateValue();
  BoundaryConditions();

  // PV to LIA y
  AdvectionY(nz_ + 1, nx_ + 1, pv_, lia_y_, &pv_new_, &lia_y_new_);
  // LIA z to SIA x
  AdvectionY(nz_, nx_ + 1, lia_z_, sia_x_, &lia_z_new_, &sia_x_new_);
  // LIA x to SIA z
  AdvectionY(nz_ + 1, nx_, lia_x_, sia_z_, &lia_x_new_, &sia_z_new_);
  // SIA y to VIA
  AdvectionY(nz_, nx_, sia_y_, via_, &sia_y_new_, &via_new_);
  UpdateValue();
  BoundaryConditions();

  // PV to LIA z
  AdvectionZ(nx_ + 1, ny_ + 1, pv_, lia_z_, &pv_new_, &lia_z_new_);
  // LIA x to SIA y
  AdvectionZ(nx_, ny_ + 1, lia_x_, sia_y_, &lia_x_new_, &sia_y_new_);
  // LIA y to SIA x
  AdvectionZ(nx_ + 1, ny_, lia_y_, sia_x_, &lia_y_new_, &sia_x_new_);
  // SIA z to VIA
  AdvectionZ(nx_, ny_, sia_z_, via_, &sia_z_new_, &via_new_);
  UpdateValue();
  BoundaryConditions();
}

#include "AdvectionX.h"
#include "AdvectionY.h"
#include "AdvectionZ.h"
