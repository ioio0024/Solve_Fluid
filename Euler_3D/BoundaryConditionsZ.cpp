#include "global_definitions.h"
#include "Fluid3D.h"

void Fluid3D::BoundaryConditionsZ() noexcept
{
  // Z Left
  // PV
  for(auto k = 0; k < GHOST + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_pv_(i, j, k) = density_pv_(i, j, GHOST);
        vx_pv_(i, j, k) = vx_pv_(i, j, GHOST);
        vy_pv_(i, j, k) = vy_pv_(i, j, GHOST);
        vz_pv_(i, j, k) = 0.0;
        pressure_pv_(i, j, k) = pressure_pv_(i, j, GHOST);
      }
    }
  }
  // LIA x
  for(auto k = 0; k < GHOST + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_lia_x_(i, j, k) = density_lia_x_(i, j, GHOST);
        vx_lia_x_(i, j, k) = vx_lia_x_(i, j, GHOST);
        vy_lia_x_(i, j, k) = vy_lia_x_(i, j, GHOST);
        vz_lia_x_(i, j, k) = 0.0;
        pressure_lia_x_(i, j, k) = pressure_lia_x_(i, j, GHOST);
      }
    }
  }
  // LIA y
  for(auto k = 0; k < GHOST + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_lia_y_(i, j, k) = density_lia_y_(i, j, GHOST);
        vx_lia_y_(i, j, k) = vx_lia_y_(i, j, GHOST);
        vy_lia_y_(i, j, k) = vy_lia_y_(i, j, GHOST);
        vz_lia_y_(i, j, k) = 0.0;
        pressure_lia_y_(i, j, k) = pressure_lia_y_(i, j, GHOST);
      }
    }
  }
  // LIA z
  for(auto k = 0; k < GHOST; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_lia_z_(i, j, k) = density_lia_z_(i, j, GHOST);
        vx_lia_z_(i, j, k) = vx_lia_z_(i, j, GHOST);
        vy_lia_z_(i, j, k) = vy_lia_z_(i, j, GHOST);
        vz_lia_z_(i, j, k) = -vz_lia_z_(i, j, GHOST);
        pressure_lia_z_(i, j, k) = pressure_lia_z_(i, j, GHOST);
      }
    }
  }
  // SIA x
  for(auto k = 0; k < GHOST; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        density_sia_x_(i, j, k) = density_sia_x_(i, j, GHOST);
        vx_sia_x_(i, j, k) = vx_sia_x_(i, j, GHOST);
        vy_sia_x_(i, j, k) = vy_sia_x_(i, j, GHOST);
        vz_sia_x_(i, j, k) = -vz_sia_x_(i, j, GHOST);
        pressure_sia_x_(i, j, k) = pressure_sia_x_(i, j, GHOST);
      }
    }
  }
  // SIA y
  for(auto k = 0; k < GHOST; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_sia_y_(i, j, k) = density_sia_y_(i, j, GHOST);
        vx_sia_y_(i, j, k) = vx_sia_y_(i, j, GHOST);
        vy_sia_y_(i, j, k) = vy_sia_y_(i, j, GHOST);
        vz_sia_y_(i, j, k) = -vz_sia_y_(i, j, GHOST);
        pressure_sia_y_(i, j, k) = pressure_sia_y_(i, j, GHOST);
      }
    }
  }
  // SIA z
  for(auto k = 0; k < GHOST + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_sia_z_(i, j, k) = density_sia_z_(i, j, GHOST);
        vx_sia_z_(i, j, k) = vx_sia_z_(i, j, GHOST);
        vy_sia_z_(i, j, k) = vy_sia_z_(i, j, GHOST);
        vz_sia_z_(i, j, k) = 0.0;
        pressure_sia_z_(i, j, k) = pressure_sia_z_(i, j, GHOST);
      }
    }
  }
  // VIA
  for(auto k = 0; k < GHOST; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_; ++i){
        density_via_(i, j, k) = density_via_(i, j, GHOST);
        vx_via_(i, j, k) = vx_via_(i, j, GHOST);
        vy_via_(i, j, k) = vy_via_(i, j, GHOST);
        vz_via_(i, j, k) = -vz_via_(i, j, GHOST);
        pressure_via_(i, j, k) = pressure_via_(i, j, GHOST);
      }
    }
  }

  // Z Right
  // PV
  for(auto k = 0; k < GHOST + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        auto n = nz_ - k;
        auto ne = nz_ - GHOST;
        density_pv_(i, j, n) = density_pv_(i, j, ne);
        vx_pv_(i, j, n) = vx_pv_(i, j, ne);
        vy_pv_(i, j, n) = vy_pv_(i, j, ne);
        vz_pv_(i, j, n) = 0.0;
        pressure_pv_(i, j, n) = pressure_pv_(i, j, ne);
      }
    }
  }
  // LIA x
  for(auto k = 0; k < GHOST + 1; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        auto n = nz_ - k;
        auto ne = nz_ - GHOST;
        density_lia_x_(i, j, n) = density_lia_x_(i, j, ne);
        vx_lia_x_(i, j, n) = vx_lia_x_(i, j, ne);
        vy_lia_x_(i, j, n) = vy_lia_x_(i, j, ne);
        vz_lia_x_(i, j, n) = 0.0;
        pressure_lia_x_(i, j, n) = pressure_lia_x_(i, j, ne);
      }
    }
  }
  // LIA y
  for(auto k = 0; k < GHOST + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        auto n = nz_ - k;
        auto ne = nz_ - GHOST;
        density_lia_y_(i, j, n) = density_lia_y_(i, j, ne);
        vx_lia_y_(i, j, n) = vx_lia_y_(i, j, ne);
        vy_lia_y_(i, j, n) = vy_lia_y_(i, j, ne);
        vz_lia_y_(i, j, n) = 0.0;
        pressure_lia_y_(i, j, n) = pressure_lia_y_(i, j, ne);
      }
    }
  }
  // LIA z
  for(auto k = 0; k < GHOST; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        auto n = nz_ - k - 1;
        auto ne = nz_ - GHOST - 1;
        density_lia_z_(i, j, n) = density_lia_z_(i, j, ne);
        vx_lia_z_(i, j, n) = vx_lia_z_(i, j, ne);
        vy_lia_z_(i, j, n) = vy_lia_z_(i, j, ne);
        vz_lia_z_(i, j, n) = -vz_lia_z_(i, j, ne);
        pressure_lia_z_(i, j, n) = pressure_lia_z_(i, j, ne);
      }
    }
  }
  // SIA x
  for(auto k = 0; k < GHOST; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_ + 1; ++i){
        auto n = nz_ - k - 1;
        auto ne = nz_ - GHOST - 1;
        density_sia_x_(i, j, n) = density_sia_x_(i, j, ne);
        vx_sia_x_(i, j, n) = vx_sia_x_(i, j, ne);
        vy_sia_x_(i, j, n) = vy_sia_x_(i, j, ne);
        vz_sia_x_(i, j, n) = -vz_sia_x_(i, j, ne);
        pressure_sia_x_(i, j, n) = pressure_sia_x_(i, j, ne);
      }
    }
  }
  // SIA y
  for(auto k = 0; k < GHOST; ++k){
    for(auto j = 0; j < ny_ + 1; ++j){
      for(auto i = 0; i < nx_; ++i){
        auto n = nz_ - k - 1;
        auto ne = nz_ - GHOST - 1;
        density_sia_y_(i, j, n) = density_sia_y_(i, j, ne);
        vx_sia_y_(i, j, n) = vx_sia_y_(i, j, ne);
        vy_sia_y_(i, j, n) = vy_sia_y_(i, j, ne);
        vz_sia_y_(i, j, n) = -vz_sia_y_(i, j, ne);
        pressure_sia_y_(i, j, n) = pressure_sia_y_(i, j, ne);
      }
    }
  }
  // SIA z
  for(auto k = 0; k < GHOST + 1; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_; ++i){
        auto n = nz_ - k;
        auto ne = nz_ - GHOST;
        density_sia_z_(i, j, n) = density_sia_z_(i, j, ne);
        vx_sia_z_(i, j, n) = vx_sia_z_(i, j, ne);
        vy_sia_z_(i, j, n) = vy_sia_z_(i, j, ne);
        vz_sia_z_(i, j, n) = 0.0;
        pressure_sia_z_(i, j, n) = pressure_sia_z_(i, j, ne);
      }
    }
  }
  // VIA
  for(auto k = 0; k < GHOST; ++k){
    for(auto j = 0; j < ny_; ++j){
      for(auto i = 0; i < nx_; ++i){
        auto n = nz_ - k - 1;
        auto ne = nz_ - GHOST - 1;
        density_via_(i, j, n) = density_via_(i, j, ne);
        vx_via_(i, j, n) = vx_via_(i, j, ne);
        vy_via_(i, j, n) = vy_via_(i, j, ne);
        vz_via_(i, j, n) = -vz_via_(i, j, ne);
        pressure_via_(i, j, n) = pressure_via_(i, j, ne);
      }
    }
  }
}
