#include "global_definitions.h"
#include "Fluid3D.h"

void Fluid3D::BoundaryConditionsZ() noexcept
{
  // Z Left
  // PV and LIA z
  bc_z_left_(0, nx_ + 1, 0, ny_ + 1, 0, GHOST + 1,
             &density_pv_,
             &vx_pv_, &vy_pv_, &vz_pv_,
             &pressure_pv_,
             &density_lia_z_,
             &vx_lia_z_, &vy_lia_z_, &vz_lia_z_,
             &pressure_lia_z_);
  // LIA x and SIA y
  bc_z_left_(0, nx_, 0, ny_ + 1, 0, GHOST + 1,
             &density_lia_x_,
             &vx_lia_x_, &vy_lia_x_, &vz_lia_x_,
             &pressure_lia_x_,
             &density_sia_y_,
             &vx_sia_y_, &vy_sia_y_, &vz_sia_y_,
             &pressure_sia_y_);
  // LIA y and SIA x
  bc_z_left_(0, nx_ + 1, 0, ny_, 0, GHOST + 1,
             &density_lia_y_,
             &vx_lia_y_, &vy_lia_y_, &vz_lia_y_,
             &pressure_lia_y_,
             &density_sia_x_,
             &vx_sia_x_, &vy_sia_x_, &vz_sia_x_,
             &pressure_sia_x_);
  // SIA z and VIA
  bc_z_left_(0, nx_, 0, ny_, 0, GHOST + 1,
             &density_sia_z_,
             &vx_sia_z_, &vy_sia_z_, &vz_sia_z_,
             &pressure_sia_z_,
             &density_via_,
             &vx_via_, &vy_via_, &vz_via_,
             &pressure_via_);

  // Z Right
  bc_z_right_(0, nx_ + 1, 0, ny_ + 1, nz_ - GHOST, nz_ + 1,
             &density_pv_,
             &vx_pv_, &vy_pv_, &vz_pv_,
             &pressure_pv_,
             &density_lia_z_,
             &vx_lia_z_, &vy_lia_z_, &vz_lia_z_,
             &pressure_lia_z_);
  // LIA x and SIA y
  bc_z_right_(0, nx_, 0, ny_ + 1, nz_ - GHOST, nz_ + 1,
             &density_lia_x_,
             &vx_lia_x_, &vy_lia_x_, &vz_lia_x_,
             &pressure_lia_x_,
             &density_sia_y_,
             &vx_sia_y_, &vy_sia_y_, &vz_sia_y_,
             &pressure_sia_y_);
  // LIA y and SIA x
  bc_z_right_(0, nx_ + 1, 0, ny_, nz_ - GHOST, nz_ + 1,
             &density_lia_y_,
             &vx_lia_y_, &vy_lia_y_, &vz_lia_y_,
             &pressure_lia_y_,
             &density_sia_x_,
             &vx_sia_x_, &vy_sia_x_, &vz_sia_x_,
             &pressure_sia_x_);
  // SIA z and VIA
  bc_z_right_(0, nx_, 0, ny_, nz_ - GHOST, nz_ + 1,
             &density_sia_z_,
             &vx_sia_z_, &vy_sia_z_, &vz_sia_z_,
             &pressure_sia_z_,
             &density_via_,
             &vx_via_, &vy_via_, &vz_via_,
             &pressure_via_);
}
