template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>
::BoundaryConditionsY() noexcept
{
  // Y Left
  // PV and LIA y
  bc_y_left_(0, nx_ + 1, 0, GHOST + 1, 0, nz_ + 1, &pv_, &lia_y_);
  // LIA z and SIA x
  bc_y_left_(0, nx_ + 1, 0, GHOST + 1, 0, nz_, &lia_z_, &sia_x_);
  // LIA x and SIA z
  bc_y_left_(0, nx_, 0, GHOST + 1, 0, nz_ + 1, &lia_x_, &sia_z_);
  // SIA y and VIA
  bc_y_left_(0, nx_, 0, GHOST + 1, 0, nz_, &sia_y_, &via_);

  // Y Right
  // PV and LIA y
  bc_y_right_(0, nx_ + 1, ny_ - GHOST, ny_ + 1, 0, nz_ + 1, &pv_, &lia_y_);
  // LIA z and SIA x
  bc_y_right_(0, nx_ + 1, ny_ - GHOST, ny_ + 1, 0, nz_, &lia_z_, &sia_x_);
  // LIA x and SIA z
  bc_y_right_(0, nx_, ny_ - GHOST, ny_ + 1, 0, nz_ + 1, &lia_x_, &sia_z_);
  // SIA y and VIA
  bc_y_right_(0, nx_, ny_ - GHOST, ny_ + 1, 0, nz_, &sia_y_, &via_);
}
