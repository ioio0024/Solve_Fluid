#include "global_definitions.h"
#include "Fluid3D.h"

void Fluid3D::BoundaryConditionsX() noexcept
{
  // X Left
  // PV
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < GHOST + 1; ++i){
        density_pv_(i, j, k) = density_pv_(GHOST, j, k);
        vx_pv_(i, j, k) = 0.0;
        vy_pv_(i, j, k) = vy_pv_(GHOST, j, k);
        vz_pv_(i, j, k) = vz_pv_(GHOST, j, k);
        pressure_pv_(i, j, k) = pressure_pv_(GHOST, j, k);
      }
    }
  }
  // LIA x
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < GHOST; ++i){
        density_lia_x_(i, j, k) = density_lia_x_(GHOST, j, k);
        vx_lia_x_(i, j, k) = -vx_lia_x_(GHOST, j, k);
        vy_lia_x_(i, j, k) = vy_lia_x_(GHOST, j, k);
        vz_lia_x_(i, j, k) = vz_lia_x_(GHOST, j, k);
        pressure_lia_x_(i, j, k) = pressure_lia_x_(GHOST, j, k);
      }
    }
  }
  // LIA y
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < GHOST + 1; ++i){
        density_lia_y_(i, j, k) = density_lia_y_(GHOST, j, k);
        vx_lia_y_(i, j, k) = 0.0;
        vy_lia_y_(i, j, k) = vy_lia_y_(GHOST, j, k);
        vz_lia_y_(i, j, k) = vz_lia_y_(GHOST, j, k);
        pressure_lia_y_(i, j, k) = pressure_lia_y_(GHOST, j, k);
      }
    }
  }
  // LIA z
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < GHOST + 1; ++i){
        density_lia_z_(i, j, k) = density_lia_z_(GHOST, j, k);
        vx_lia_z_(i, j, k) = 0.0;
        vy_lia_z_(i, j, k) = vy_lia_z_(GHOST, j, k);
        vz_lia_z_(i, j, k) = vz_lia_z_(GHOST, j, k);
        pressure_lia_z_(i, j, k) = pressure_lia_z_(GHOST, j, k);
      }
    }
  }
  // SIA x
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < GHOST + 1; ++i){
        density_sia_x_(i, j, k) = density_sia_x_(GHOST, j, k);
        vx_sia_x_(i, j, k) = 0.0;
        vy_sia_x_(i, j, k) = vy_sia_x_(GHOST, j, k);
        vz_sia_x_(i, j, k) = vz_sia_x_(GHOST, j, k);
        pressure_sia_x_(i, j, k) = pressure_sia_x_(GHOST, j, k);
      }
    }
  }
  // SIA y
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < GHOST; ++i){
        density_sia_y_(i, j, k) = density_sia_y_(GHOST, j, k);
        vx_sia_y_(i, j, k) = -vx_sia_y_(GHOST, j, k);
        vy_sia_y_(i, j, k) = vy_sia_y_(GHOST, j, k);
        vz_sia_y_(i, j, k) = vz_sia_y_(GHOST, j, k);
        pressure_sia_y_(i, j, k) = pressure_sia_y_(GHOST, j, k);
      }
    }
  }
  // SIA z
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < GHOST; ++i){
        density_sia_z_(i, j, k) = density_sia_z_(GHOST, j, k);
        vx_sia_z_(i, j, k) = -vx_sia_z_(GHOST, j, k);
        vy_sia_z_(i, j, k) = vy_sia_z_(GHOST, j, k);
        vz_sia_z_(i, j, k) = vz_sia_z_(GHOST, j, k);
        pressure_sia_z_(i, j, k) = pressure_sia_z_(GHOST, j, k);
      }
    }
  }
  // VIA
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < GHOST; ++i){
        density_via_(i, j, k) = density_via_(GHOST, j, k);
        vx_via_(i, j, k) = -vx_via_(GHOST, j, k);
        vy_via_(i, j, k) = vy_via_(GHOST, j, k);
        vz_via_(i, j, k) = vz_via_(GHOST, j, k);
        pressure_via_(i, j, k) = pressure_via_(GHOST, j, k);
      }
    }
  }

  // X Right
  // PV
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < GHOST + 1; ++i){
        auto n = nx_ - i;
        auto ne = nx_ - GHOST;
        density_pv_(n, j, k) = density_pv_(ne, j, k);
        vx_pv_(n, j, k) = 0.0;
        vy_pv_(n, j, k) = vy_pv_(ne, j, k);
        vz_pv_(n, j, k) = vz_pv_(ne, j, k);
        pressure_pv_(n, j, k) = pressure_pv_(ne, j, k);
      }
    }
  }
  // LIA x
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < GHOST; ++i){
        auto n = nx_ - i - 1;
        auto ne = nx_ - GHOST - 1;
        density_lia_x_(n, j, k) = density_lia_x_(ne, j, k);
        vx_lia_x_(n, j, k) = -vx_lia_x_(ne, j, k);
        vy_lia_x_(n, j, k) = vy_lia_x_(ne, j, k);
        vz_lia_x_(n, j, k) = vz_lia_x_(ne, j, k);
        pressure_lia_x_(n, j, k) = pressure_lia_x_(ne, j, k);
      }
    }
  }
  // LIA y
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < GHOST + 1; ++i){
        auto n = nx_ - i;
        auto ne = nx_ - GHOST;
        density_lia_y_(n, j, k) = density_lia_y_(ne, j, k);
        vx_lia_y_(n, j, k) = 0.0;
        vy_lia_y_(n, j, k) = vy_lia_y_(ne, j, k);
        vz_lia_y_(n, j, k) = vz_lia_y_(ne, j, k);
        pressure_lia_y_(n, j, k) = pressure_lia_y_(ne, j, k);
      }
    }
  }
  // LIA z
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < GHOST + 1; ++i){
        auto n = nx_ - i;
        auto ne = nx_ - GHOST;
        density_lia_z_(n, j, k) = density_lia_z_(ne, j, k);
        vx_lia_z_(n, j, k) = 0.0;
        vy_lia_z_(n, j, k) = vy_lia_z_(ne, j, k);
        vz_lia_z_(n, j, k) = vz_lia_z_(ne, j, k);
        pressure_lia_z_(n, j, k) = pressure_lia_z_(ne, j, k);
      }
    }
  }
  // SIA x
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < GHOST + 1; ++i){
        auto n = nx_ - i;
        auto ne = nx_ - GHOST;
        density_sia_x_(n, j, k) = density_sia_x_(ne, j, k);
        vx_sia_x_(n, j, k) = 0.0;
        vy_sia_x_(n, j, k) = vy_sia_x_(ne, j, k);
        vz_sia_x_(n, j, k) = vz_sia_x_(ne, j, k);
        pressure_sia_x_(n, j, k) = pressure_sia_x_(ne, j, k);
      }
    }
  }
  // SIA y
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < GHOST; ++i){
        auto n = nx_ - i - 1;
        auto ne = nx_ - GHOST - 1;
        density_sia_y_(n, j, k) = density_sia_y_(ne, j, k);
        vx_sia_y_(n, j, k) = -vx_sia_y_(ne, j, k);
        vy_sia_y_(n, j, k) = vy_sia_y_(ne, j, k);
        vz_sia_y_(n, j, k) = vz_sia_y_(ne, j, k);
        pressure_sia_y_(n, j, k) = pressure_sia_y_(ne, j, k);
      }
    }
  }
  // SIA z
  for(auto k = 0; k < nz_ + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < GHOST; ++i){
        auto n = nx_ - i - 1;
        auto ne = nx_ - GHOST - 1;
        density_sia_z_(n, j, k) = density_sia_z_(ne, j, k);
        vx_sia_z_(n, j, k) = -vx_sia_z_(ne, j, k);
        vy_sia_z_(n, j, k) = vy_sia_z_(ne, j, k);
        vz_sia_z_(n, j, k) = vz_sia_z_(ne, j, k);
        pressure_sia_z_(n, j, k) = pressure_sia_z_(ne, j, k);
      }
    }
  }
  // VIA
  for(auto k = 0; k < nz_; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < GHOST; ++i){
        auto n = nx_ - i - 1;
        auto ne = nx_ - GHOST - 1;
        density_via_(n, j, k) = density_via_(ne, j, k);
        vx_via_(n, j, k) = -vx_via_(ne, j, k);
        vy_via_(n, j, k) = vy_via_(ne, j, k);
        vz_via_(n, j, k) = vz_via_(ne, j, k);
        pressure_via_(n, j, k) = pressure_via_(ne, j, k);
      }
    }
  }
}
