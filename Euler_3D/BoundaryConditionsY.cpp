#include "global_definitions.h"
#include "Fluid3D.h"

void Fluid3D::BoundaryConditionsY() noexcept
{
  // Y Left
  // PV
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < GHOST + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_pv_(i, j, k) = density_pv_(i, GHOST, k);
        vx_pv_(i, j, k) = vx_pv_(i, GHOST, k);
        vy_pv_(i, j, k) = 0.0;
        vz_pv_(i, j, k) = vz_pv_(i, GHOST, k);
        pressure_pv_(i, j, k) = pressure_pv_(i, GHOST, k);
      }
    }
  }
  // LIA x
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < GHOST + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_lia_x_(i, j, k) = density_lia_x_(i, GHOST, k);
        vx_lia_x_(i, j, k) = vx_lia_x_(i, GHOST, k);
        vy_lia_x_(i, j, k) = 0.0;
        vz_lia_x_(i, j, k) = vz_lia_x_(i, GHOST, k);
        pressure_lia_x_(i, j, k) = pressure_lia_x_(i, GHOST, k);
      }
    }
  }
  // LIA y
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < GHOST; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_lia_y_(i, j, k) = density_lia_y_(i, GHOST, k);
        vx_lia_y_(i, j, k) = vx_lia_y_(i, GHOST, k);
        vy_lia_y_(i, j, k) = -vy_lia_y_(i, GHOST, k);
        vz_lia_y_(i, j, k) = vz_lia_y_(i, GHOST, k);
        pressure_lia_y_(i, j, k) = pressure_lia_y_(i, GHOST, k);
      }
    }
  }
  // LIA z
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < GHOST + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_lia_z_(i, j, k) = density_lia_z_(i, GHOST, k);
        vx_lia_z_(i, j, k) = vx_lia_z_(i, GHOST, k);
        vy_lia_z_(i, j, k) = 0.0;
        vz_lia_z_(i, j, k) = vz_lia_z_(i, GHOST, k);
        pressure_lia_z_(i, j, k) = pressure_lia_z_(i, GHOST, k);
      }
    }
  }
  // SIA x
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < GHOST; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_sia_x_(i, j, k) = density_sia_x_(i, GHOST, k);
        vx_sia_x_(i, j, k) = vx_sia_x_(i, GHOST, k);
        vy_sia_x_(i, j, k) = -vy_sia_x_(i, GHOST, k);
        vz_sia_x_(i, j, k) = vz_sia_x_(i, GHOST, k);
        pressure_sia_x_(i, j, k) = pressure_sia_x_(i, GHOST, k);
      }
    }
  }
  // SIA y
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < GHOST + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_sia_y_(i, j, k) = density_sia_y_(i, GHOST, k);
        vx_sia_y_(i, j, k) = vx_sia_y_(i, GHOST, k);
        vy_sia_y_(i, j, k) = 0.0;
        vz_sia_y_(i, j, k) = vz_sia_y_(i, GHOST, k);
        pressure_sia_y_(i, j, k) = pressure_sia_y_(i, GHOST, k);
      }
    }
  }
  // SIA z
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < GHOST; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_sia_z_(i, j, k) = density_sia_z_(i, GHOST, k);
        vx_sia_z_(i, j, k) = vx_sia_z_(i, GHOST, k);
        vy_sia_z_(i, j, k) = -vy_sia_z_(i, GHOST, k);
        vz_sia_z_(i, j, k) = vz_sia_z_(i, GHOST, k);
        pressure_sia_z_(i, j, k) = pressure_sia_z_(i, GHOST, k);
      }
    }
  }
  // VIA
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < GHOST; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_via_(i, j, k) = density_via_(i, GHOST, k);
        vx_via_(i, j, k) = vx_via_(i, GHOST, k);
        vy_via_(i, j, k) = -vy_via_(i, GHOST, k);
        vz_via_(i, j, k) = vz_via_(i, GHOST, k);
        pressure_via_(i, j, k) = pressure_via_(i, GHOST, k);
      }
    }
  }

  // Y Right
  // PV
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < GHOST + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        auto n = ny_ - j;
        auto ne = ny_ - GHOST;
        density_pv_(i, n, k) = density_pv_(i, ne, k);
        vx_pv_(i, n, k) = vx_pv_(i, ne, k);
        vy_pv_(i, n, k) = 0.0;
        vz_pv_(i, n, k) = vz_pv_(i, ne, k);
        pressure_pv_(i, n, k) = pressure_pv_(i, ne, k);
      }
    }
  }
  // LIA x
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < GHOST + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        auto n = ny_ - j;
        auto ne = ny_ - GHOST;
        density_lia_x_(i, n, k) = density_lia_x_(i, ne, k);
        vx_lia_x_(i, n, k) = vx_lia_x_(i, ne, k);
        vy_lia_x_(i, n, k) = 0.0;
        vz_lia_x_(i, n, k) = vz_lia_x_(i, ne, k);
        pressure_lia_x_(i, n, k) = pressure_lia_x_(i, ne, k);
      }
    }
  }
  // LIA y
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < GHOST; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        auto n = ny_ - j - 1;
        auto ne = ny_ - GHOST - 1;
        density_lia_y_(i, n, k) = density_lia_y_(i, ne, k);
        vx_lia_y_(i, n, k) = vx_lia_y_(i, ne, k);
        vy_lia_y_(i, n, k) = -vy_lia_y_(i, ne, k);
        vz_lia_y_(i, n, k) = vz_lia_y_(i, ne, k);
        pressure_lia_y_(i, n, k) = pressure_lia_y_(i, ne, k);
      }
    }
  }
  // LIA z
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < GHOST + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        auto n = ny_ - j;
        auto ne = ny_ - GHOST;
        density_lia_z_(i, n, k) = density_lia_z_(i, ne, k);
        vx_lia_z_(i, n, k) = vx_lia_z_(i, ne, k);
        vy_lia_z_(i, n, k) = 0.0;
        vz_lia_z_(i, n, k) = vz_lia_z_(i, ne, k);
        pressure_lia_z_(i, n, k) = pressure_lia_z_(i, ne, k);
      }
    }
  }
  // SIA x
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < GHOST; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        auto n = ny_ - j - 1;
        auto ne = ny_ - GHOST - 1;
        density_sia_x_(i, n, k) = density_sia_x_(i, ne, k);
        vx_sia_x_(i, n, k) = vx_sia_x_(i, ne, k);
        vy_sia_x_(i, n, k) = -vy_sia_x_(i, ne, k);
        vz_sia_x_(i, n, k) = vz_sia_x_(i, ne, k);
        pressure_sia_x_(i, n, k) = pressure_sia_x_(i, ne, k);
      }
    }
  }
  // SIA y
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < GHOST + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        auto n = ny_ - j;
        auto ne = ny_ - GHOST;
        density_sia_y_(i, n, k) = density_sia_y_(i, ne, k);
        vx_sia_y_(i, n, k) = vx_sia_y_(i, ne, k);
        vy_sia_y_(i, n, k) = 0.0;
        vz_sia_y_(i, n, k) = vz_sia_y_(i, ne, k);
        pressure_sia_y_(i, n, k) = pressure_sia_y_(i, ne, k);
      }
    }
  }
  // SIA z
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < GHOST; ++j){
      for(auto i = 0; i < nx_; ++i){
        auto n = ny_ - j - 1;
        auto ne = ny_ - GHOST - 1;
        density_sia_z_(i, n, k) = density_sia_z_(i, ne, k);
        vx_sia_z_(i, n, k) = vx_sia_z_(i, ne, k);
        vy_sia_z_(i, n, k) = -vy_sia_z_(i, ne, k);
        vz_sia_z_(i, n, k) = vz_sia_z_(i, ne, k);
        pressure_sia_z_(i, n, k) = pressure_sia_z_(i, ne, k);
      }
    }
  }
  // VIA
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < GHOST; ++j){
      for(auto i = 0; i < nx_; ++i){
        auto n = ny_ - j - 1;
        auto ne = ny_ - GHOST - 1;
        density_via_(i, n, k) = density_via_(i, ne, k);
        vx_via_(i, n, k) = vx_via_(i, ne, k);
        vy_via_(i, n, k) = -vy_via_(i, ne, k);
        vz_via_(i, n, k) = vz_via_(i, ne, k);
        pressure_via_(i, n, k) = pressure_via_(i, ne, k);
      }
    }
  }
}
