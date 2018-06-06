#include "global_definitions.h"
#include "Fluid3D.h"

void Fluid3D::Advection()
{
  // PV to LIA x
  AdvectionX(ny_ + 1, nz_ + 1,
             density_pv_, vx_pv_, vy_pv_,
             vz_pv_, pressure_pv_,
             density_lia_x_, vx_lia_x_, vy_lia_x_,
             vz_lia_x_, pressure_lia_x_,
             &density_pv_new_, &vx_pv_new_, &vy_pv_new_,
             &vz_pv_new_, &pressure_pv_new_,
             &density_lia_x_new_, &vx_lia_x_new_, &vy_lia_x_new_,
             &vz_lia_x_new_, &pressure_lia_x_new_);
  // LIA y to SIA z
  AdvectionX(ny_, nz_ + 1,
             density_lia_y_, vx_lia_y_, vy_lia_y_,
             vz_lia_y_, pressure_lia_y_,
             density_sia_z_, vx_sia_z_, vy_sia_z_,
             vz_sia_z_, pressure_sia_z_,
             &density_lia_y_new_, &vx_lia_y_new_, &vy_lia_y_new_,
             &vz_lia_y_new_, &pressure_lia_y_new_,
             &density_sia_z_new_, &vx_sia_z_new_, &vy_sia_z_new_,
             &vz_sia_z_new_, &pressure_sia_z_new_);
  // LIA z to SIA y
  AdvectionX(ny_ + 1, nz_,
             density_lia_z_, vx_lia_z_, vy_lia_z_,
             vz_lia_z_, pressure_lia_z_,
             density_sia_y_, vx_sia_y_, vy_sia_y_,
             vz_sia_y_, pressure_sia_y_,
             &density_lia_z_new_, &vx_lia_z_new_, &vy_lia_z_new_,
             &vz_lia_z_new_, &pressure_lia_z_new_,
             &density_sia_y_new_, &vx_sia_y_new_, &vy_sia_y_new_,
             &vz_sia_y_new_, &pressure_sia_y_new_);
  // SIA x to VIA
  AdvectionX(ny_, nz_,
             density_sia_x_, vx_sia_x_, vy_sia_x_,
             vz_sia_x_, pressure_sia_x_,
             density_via_, vx_via_, vy_via_,
             vz_via_, pressure_via_,
             &density_sia_x_new_, &vx_sia_x_new_, &vy_sia_x_new_,
             &vz_sia_x_new_, &pressure_sia_x_new_,
             &density_via_new_, &vx_via_new_, &vy_via_new_,
             &vz_via_new_, &pressure_via_new_);
  UpdateValue();
  BoundaryConditions();

  // PV to LIA y
  AdvectionY(nz_ + 1, nx_ + 1,
             density_pv_, vx_pv_, vy_pv_,
             vz_pv_, pressure_pv_,
             density_lia_y_, vx_lia_y_, vy_lia_y_,
             vz_lia_y_, pressure_lia_y_,
             &density_pv_new_, &vx_pv_new_, &vy_pv_new_,
             &vz_pv_new_, &pressure_pv_new_,
             &density_lia_y_new_, &vx_lia_y_new_, &vy_lia_y_new_,
             &vz_lia_y_new_, &pressure_lia_y_new_);
  // LIA z to SIA x
  AdvectionY(nz_, nx_ + 1,
             density_lia_z_, vx_lia_z_, vy_lia_z_,
             vz_lia_z_, pressure_lia_z_,
             density_sia_x_, vx_sia_x_, vy_sia_x_,
             vz_sia_x_, pressure_sia_x_,
             &density_lia_z_new_, &vx_lia_z_new_, &vy_lia_z_new_,
             &vz_lia_z_new_,&pressure_lia_z_new_,
             &density_sia_x_new_, &vx_sia_x_new_, &vy_sia_x_new_,
             &vz_sia_x_new_, &pressure_sia_x_new_);
  // LIA x to SIA z
  AdvectionY(nz_ + 1, nx_,
             density_lia_x_, vx_lia_x_, vy_lia_x_,
             vz_lia_x_, pressure_lia_x_,
             density_sia_z_, vx_sia_z_, vy_sia_z_,
             vz_sia_z_, pressure_sia_z_,
             &density_lia_x_new_, &vx_lia_x_new_, &vy_lia_x_new_,
             &vz_lia_x_new_, &pressure_lia_x_new_,
             &density_sia_z_new_, &vx_sia_z_new_, &vy_sia_z_new_,
             &vz_sia_z_new_, &pressure_sia_z_new_);
  // SIA y to VIA
  AdvectionY(nz_, nx_,
             density_sia_y_, vx_sia_y_, vy_sia_y_,
             vz_sia_y_, pressure_sia_y_,
             density_via_, vx_via_, vy_via_,
             vz_via_, pressure_via_,
             &density_sia_y_new_, &vx_sia_y_new_, &vy_sia_y_new_,
             &vz_sia_y_new_, &pressure_sia_y_new_,
             &density_via_new_, &vx_via_new_, &vy_via_new_,
             &vz_via_new_, &pressure_via_new_);
  UpdateValue();
  BoundaryConditions();

  // PV to LIA z
  AdvectionZ(nx_ + 1, ny_ + 1,
             density_pv_, vx_pv_, vy_pv_,
             vz_pv_, pressure_pv_,
             density_lia_z_, vx_lia_z_, vy_lia_z_,
             vz_lia_z_, pressure_lia_z_,
             &density_pv_new_, &vx_pv_new_, &vy_pv_new_,
             &vz_pv_new_, &pressure_pv_new_,
             &density_lia_z_new_, &vx_lia_z_new_, &vy_lia_z_new_,
             &vz_lia_z_new_, &pressure_lia_z_new_);
  // LIA x to SIA y
  AdvectionZ(nx_, ny_ + 1,
             density_lia_x_, vx_lia_x_, vy_lia_x_,
             vz_lia_x_, pressure_lia_x_,
             density_sia_y_, vx_sia_y_, vy_sia_y_,
             vz_sia_y_, pressure_sia_y_,
             &density_lia_x_new_, &vx_lia_x_new_, &vy_lia_x_new_,
             &vz_lia_x_new_, &pressure_lia_x_new_,
             &density_sia_y_new_, &vx_sia_y_new_, &vy_sia_y_new_,
             &vz_sia_y_new_, &pressure_sia_y_new_);
  // LIA y to SIA x
  AdvectionZ(nx_ + 1, ny_,
             density_lia_y_, vx_lia_y_, vy_lia_y_,
             vz_lia_y_, pressure_lia_y_,
             density_sia_x_, vx_sia_x_, vy_sia_x_,
             vz_sia_x_, pressure_sia_x_,
             &density_lia_y_new_, &vx_lia_y_new_, &vy_lia_y_new_,
             &vz_lia_y_new_, &pressure_lia_y_new_,
             &density_sia_x_new_, &vx_sia_x_new_, &vy_sia_x_new_,
             &vz_sia_x_new_, &pressure_sia_x_new_);
  // SIA z to VIA
  AdvectionZ(nx_, ny_,
             density_sia_z_, vx_sia_z_, vy_sia_z_,
             vz_sia_z_, pressure_sia_z_,
             density_via_, vx_via_, vy_via_,
             vz_via_, pressure_via_,
             &density_sia_z_new_, &vx_sia_z_new_, &vy_sia_z_new_,
             &vz_sia_z_new_, &pressure_sia_z_new_,
             &density_via_new_, &vx_via_new_, &vy_via_new_,
             &vz_via_new_, &pressure_via_new_);
  UpdateValue();
  BoundaryConditions();
}
