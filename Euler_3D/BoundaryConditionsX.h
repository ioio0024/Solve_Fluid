template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
::BoundaryConditionsX() noexcept
{
  // X Left
  // PV and LIA x
  bc_x_left_(0, GHOST + 1, 0, ny_ + 1, 0, nz_ + 1, &pv_, &lia_x_);
  // LIA y and SIA z
  bc_x_left_(0, GHOST + 1, 0, ny_, 0, nz_ + 1, &lia_y_, &sia_z_);
  // LIA z snd SIA y
  bc_x_left_(0, GHOST + 1, 0, ny_ + 1, 0, nz_, &lia_z_, &sia_y_);
  // SIA x and VIA
  bc_x_left_(0, GHOST + 1, 0, ny_, 0, nz_, &sia_x_, &via_);

  // X Right
  // PV and LIA x
  bc_x_right_(nx_ - GHOST, nx_ + 1, 0, ny_ + 1, 0, nz_ + 1, &pv_, &lia_x_);
  // LIA y and SIA z
  bc_x_right_(nx_ - GHOST, nx_ + 1, 0, ny_, 0, nz_ + 1, &lia_y_, &sia_z_);
  // LIA z and SIA y
  bc_x_right_(nx_ - GHOST, nx_ + 1, 0, ny_ + 1, 0, nz_, &lia_z_, &sia_y_);
  // SIA x and VIA
  bc_x_right_(nx_ - GHOST, nx_ + 1, 0, ny_, 0, nz_, &sia_x_, &via_);
}
