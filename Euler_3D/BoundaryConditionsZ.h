template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
::BoundaryConditionsZ() noexcept
{
  // Z Left
  // PV and LIA z
  bc_z_left_(0, nx_ + 1, 0, ny_ + 1, 0, GHOST + 1, &pv_, &lia_z_);
  // LIA x and SIA y
  bc_z_left_(0, nx_, 0, ny_ + 1, 0, GHOST + 1, &lia_x_, &sia_y_);
  // LIA y and SIA x
  bc_z_left_(0, nx_ + 1, 0, ny_, 0, GHOST + 1, &lia_y_, &sia_x_);
  // SIA z and VIA
  bc_z_left_(0, nx_, 0, ny_, 0, GHOST + 1, &sia_z_, &via_);

  // Z Right
  bc_z_right_(0, nx_ + 1, 0, ny_ + 1, nz_ - GHOST, nz_ + 1, &pv_, &lia_z_);
  // LIA x and SIA y
  bc_z_right_(0, nx_, 0, ny_ + 1, nz_ - GHOST, nz_ + 1, &lia_x_, &sia_y_);
  // LIA y and SIA x
  bc_z_right_(0, nx_ + 1, 0, ny_, nz_ - GHOST, nz_ + 1, &lia_y_, &sia_x_);
  // SIA z and VIA
  bc_z_right_(0, nx_, 0, ny_, nz_ - GHOST, nz_ + 1, &sia_z_, &via_);
}
