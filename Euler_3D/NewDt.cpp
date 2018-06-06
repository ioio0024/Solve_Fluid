#include "global_definitions.h"
#include "Fluid3D.h"

#include <cmath>

void Fluid3D::NewDt() noexcept
{
  for(auto k = GHOST; k < nz_ + 1 - GHOST; ++k){
    Real vz, cs;
    auto dz = z_[k + 1] - z_[k];
    for(auto j = GHOST; j < ny_ + 1 - GHOST; ++j){
      Real vy;
      auto dy = y_[j + 1] - y_[j];
      for(auto i = GHOST; i < nx_ + 1 - GHOST; ++i){
        auto dx = x_[i + 1] - x_[i];
        auto vx = vx_pv_(i, j, k);
        vy = vy_pv_(i, j, k);
        vz = vz_pv_(i, j, k);
        cs = std::sqrt(gamma_ * pressure_pv_(i, j, k) / density_pv_(i, j, k));

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
