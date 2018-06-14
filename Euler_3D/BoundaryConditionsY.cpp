#include "global_definitions.h"
#include "Fluid3D.h"

void Fluid3D::BoundaryConditionsY() noexcept
{
  // Y Left
  // PV and LIA y
  bc_y_left_(0, nx_ + 1, 0, GHOST + 1, 0, nz_ + 1,
             &density_pv_,
             &vx_pv_, &vy_pv_, &vz_pv_,
             &pressure_pv_,
             &density_lia_y_,
             &vx_lia_y_, &vy_lia_y_, &vz_lia_y_,
             &pressure_lia_y_);
  // LIA z and SIA x
  bc_y_left_(0, nx_ + 1, 0, GHOST + 1, 0, nz_,
             &density_lia_z_,
             &vx_lia_z_, &vy_lia_z_, &vz_lia_z_,
             &pressure_lia_z_,
             &density_sia_x_,
             &vx_sia_x_, &vy_sia_x_, &vz_sia_x_,
             &pressure_sia_x_);
  // LIA x and SIA z
  bc_y_left_(0, nx_, 0, GHOST + 1, 0, nz_ + 1,
             &density_lia_x_,
             &vx_lia_x_, &vy_lia_x_, &vz_lia_x_,
             &pressure_lia_x_,
             &density_sia_z_,
             &vx_sia_z_, &vy_sia_z_, &vz_sia_z_,
             &pressure_sia_z_);
  // SIA y and VIA
  bc_y_left_(0, nx_, 0, GHOST + 1, 0, nz_,
             &density_sia_y_,
             &vx_sia_y_, &vy_sia_y_, &vz_sia_y_,
             &pressure_sia_y_,
             &density_via_,
             &vx_via_, &vy_via_, &vz_via_,
             &pressure_via_);

  // Y Right
  // PV and LIA y
  bc_y_right_(0, nx_ + 1, ny_ - GHOST, ny_ + 1, 0, nz_ + 1,
              &density_pv_,
              &vx_pv_, &vy_pv_, &vz_pv_,
              &pressure_pv_,
              &density_lia_y_,
              &vx_lia_y_, &vy_lia_y_, &vz_lia_y_,
              &pressure_lia_y_);
  // LIA z and SIA x
  bc_y_right_(0, nx_ + 1, ny_ - GHOST, ny_ + 1, 0, nz_,
              &density_lia_z_,
              &vx_lia_z_, &vy_lia_z_, &vz_lia_z_,
              &pressure_lia_z_,
              &density_sia_x_,
              &vx_sia_x_, &vy_sia_x_, &vz_sia_x_,
              &pressure_sia_x_);
  // LIA x and SIA z
  bc_y_right_(0, nx_, ny_ - GHOST, ny_ + 1, 0, nz_ + 1,
              &density_lia_x_,
              &vx_lia_x_, &vy_lia_x_, &vz_lia_x_,
              &pressure_lia_x_,
              &density_sia_z_,
              &vx_sia_z_, &vy_sia_z_, &vz_sia_z_,
              &pressure_sia_z_);
  // SIA y and VIA
  bc_y_right_(0, nx_, ny_ - GHOST, ny_ + 1, 0, nz_,
              &density_sia_y_,
              &vx_sia_y_, &vy_sia_y_, &vz_sia_y_,
              &pressure_sia_y_,
              &density_via_,
              &vx_via_, &vy_via_, &vz_via_,
              &pressure_via_);
}
