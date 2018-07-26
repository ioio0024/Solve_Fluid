template <typename BCxl, typename BCxr,
          typename BCyl, typename BCyr,
          typename BCzl, typename BCzr,
          typename CIPx, typename CIPy, typename CIPz>
void Fluid3D<BCxl, BCxr, BCyl, BCyr, BCzl, BCzr, CIPx, CIPy, CIPz>::NewDt() noexcept
{
  for(auto k = GHOST; k < nz_ + 1 - GHOST; ++k){
    Real vz, cs;
    auto dz = z_(k + 1) - z_(k);
    for(auto j = GHOST; j < ny_ + 1 - GHOST; ++j){
      Real vy;
      auto dy = y_(j + 1) - y_(j);
      for(auto i = GHOST; i < nx_ + 1 - GHOST; ++i){
        auto dx = x_(i + 1) - x_(i);
        auto vx = pv_.u(i, j, k);
        vy = pv_.v(i, j, k);
        vz = pv_.w(i, j, k);
        cs = std::sqrt(gamma_ * pv_.p(i, j, k) / pv_.d(i, j, k));

        Real vmax = 1.0e-10;
        vmax = std::fmax(vmax, std::fabs(vx));
        vmax = std::fmax(vmax, std::fabs(vx + cs));
        vmax = std::fmax(vmax, std::fabs(vx - cs));
        dt_ = std::fmin(dt_, CFL * dx / (vmax + EPS));

        vmax = 1.0e-10;
        vmax = std::fmax(vmax, std::fabs(vy));
        vmax = std::fmax(vmax, std::fabs(vy + cs));
        vmax = std::fmax(vmax, std::fabs(vy - cs));
        dt_ = std::fmin(dt_, CFL * dy / (vmax + EPS));

        vmax = 1.0e-10;
        vmax = std::fmax(vmax, std::fabs(vz));
        vmax = std::fmax(vmax, std::fabs(vz + cs));
        vmax = std::fmax(vmax, std::fabs(vz - cs));
        dt_ = std::fmin(dt_, CFL * dz / (vmax + EPS));
      }
    }
  }
}
